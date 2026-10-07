import numpy as np
import os


class XFeat:
    def __init__(
        self,
        top_k: int = 4096,
        detection_threshold: float = 0.05,
        intra_op_threads: int | None = 2,
    ) -> None:
        if top_k < 1:
            raise ValueError("top_k must be positive.")
        if not 0.0 <= detection_threshold <= 1.0:
            raise ValueError("detection_threshold must be in [0, 1].")
        if intra_op_threads is not None and intra_op_threads < 1:
            raise ValueError("intra_op_threads must be positive or None.")

        model_path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)), "models", "xfeat_backbone.onnx"
        )
        if not os.path.isfile(model_path):
            raise FileNotFoundError(
                f"XFeat ONNX model does not exist: {model_path}"
            )
        import onnxruntime as ort

        options = ort.SessionOptions()
        if intra_op_threads is not None:
            options.intra_op_num_threads = intra_op_threads
        self.session = ort.InferenceSession(
            model_path, sess_options=options, providers=["CPUExecutionProvider"]
        )
        if len(self.session.get_inputs()) != 1 or len(self.session.get_outputs()) != 3:
            raise ValueError(
                "Expected an XFeat backbone with one input and three outputs "
                "(features, keypoint logits, reliability map)."
            )
        self.input_name = self.session.get_inputs()[0].name
        self.top_k = top_k
        self.detection_threshold = detection_threshold

    def detect_and_compute(
        self,
        image: np.ndarray,
        top_k: int | None = None,
        detection_threshold: float | None = None,
    ) -> dict[str, np.ndarray]:
        top_k = self.top_k if top_k is None else top_k
        threshold = (
            self.detection_threshold
            if detection_threshold is None
            else detection_threshold
        )
        if top_k < 1:
            raise ValueError("top_k must be positive.")
        if not 0.0 <= threshold <= 1.0:
            raise ValueError("detection_threshold must be in [0, 1].")

        tensor = np.expand_dims(np.moveaxis(image, -1, 0), axis=0).astype(
            np.float32, copy=False
        )

        outputs = self.session.run(None, {self.input_name: tensor})
        descriptors, heatmap, reliability = _order_outputs(
            outputs, [output.name for output in self.session.get_outputs()]
        )
        points = nonmax_suppress(heatmap, threshold)
        if not len(points):
            return {
                "features": np.empty((0, 3), dtype=np.float32),
                "descriptors": np.empty((0, descriptors.shape[1]), dtype=np.float32),
            }

        scores = _sample_nearest(heatmap, points) * _sample_bilinear(
            reliability, points
        )
        order = np.argsort(-scores, kind="stable")[:top_k]
        points = points[order]
        descriptors = _sample_bicubic(descriptors[0], points / 8.0)
        descriptors /= np.maximum(
            np.linalg.norm(descriptors, axis=1, keepdims=True),
            np.finfo(np.float32).eps,
        )

        points = np.column_stack((points, np.ones(len(points), dtype=np.float32)))
        return {
            "features": points.astype(np.float32, copy=False),
            "descriptors": descriptors.astype(np.float32, copy=False),
        }


def _order_outputs(
    outputs: list[np.ndarray], names: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    descriptor = next((value for value in outputs if value.shape[1] == 64), None)
    maps = {
        name: value[0, 0]
        for name, value in zip(names, outputs)
        if value.shape[1] == 1
    }
    heatmap = maps.get("heatmap")
    reliability = maps.get("reliability")
    if descriptor is None or heatmap is None or reliability is None:
        raise ValueError("The ONNX model outputs do not match the XFeat backbone.")
    return descriptor, heatmap, reliability


def nonmax_suppress(heatmap: np.ndarray, threshold: float) -> np.ndarray:
    padded = np.pad(heatmap, 2, constant_values=-np.inf)
    height, width = heatmap.shape
    horizontal_maximum = padded[:, :width]
    for offset in range(1, 5):
        horizontal_maximum = np.maximum(
            horizontal_maximum, padded[:, offset : offset + width]
        )
    maximum = horizontal_maximum[:height]
    for offset in range(1, 5):
        maximum = np.maximum(
            maximum, horizontal_maximum[offset : offset + height]
        )
    y, x = np.nonzero((heatmap == maximum) & (heatmap > threshold))
    return np.column_stack((x, y)).astype(np.float32, copy=False)


def _sample_nearest(map_: np.ndarray, points: np.ndarray) -> np.ndarray:
    x = np.clip(np.rint(points[:, 0]).astype(np.intp), 0, map_.shape[1] - 1)
    y = np.clip(np.rint(points[:, 1]).astype(np.intp), 0, map_.shape[0] - 1)
    return map_[y, x]


def _sample_bilinear(map_: np.ndarray, points: np.ndarray) -> np.ndarray:
    height, width = map_.shape[-2:]
    x = np.clip(points[:, 0], 0, width - 1)
    y = np.clip(points[:, 1], 0, height - 1)
    x0, y0 = x.astype(np.intp), y.astype(np.intp)
    x1, y1 = np.minimum(x0 + 1, width - 1), np.minimum(y0 + 1, height - 1)
    wx, wy = x - x0, y - y0
    return (
        (1 - wx) * (1 - wy) * map_[..., y0, x0]
        + wx * (1 - wy) * map_[..., y0, x1]
        + (1 - wx) * wy * map_[..., y1, x0]
        + wx * wy * map_[..., y1, x1]
    )


def _sample_bicubic(map_: np.ndarray, points: np.ndarray) -> np.ndarray:
    # XFeat uses bicubic descriptor interpolation. Keys' cubic kernel with a=-0.75
    # matches PyTorch's bicubic sampling closely while requiring only NumPy.
    x, y = points[:, 0], points[:, 1]
    x_base, y_base = np.floor(x).astype(np.intp), np.floor(y).astype(np.intp)
    offsets = np.arange(-1, 3, dtype=np.intp)
    x_indices = np.clip(x_base[:, None] + offsets, 0, map_.shape[2] - 1)
    y_indices = np.clip(y_base[:, None] + offsets, 0, map_.shape[1] - 1)
    x_weights = _cubic(x[:, None] - (x_base[:, None] + offsets))
    y_weights = _cubic(y[:, None] - (y_base[:, None] + offsets))

    flat_map = map_.transpose(1, 2, 0).reshape(-1, map_.shape[0])
    flat_indices = (
        y_indices[:, :, None] * map_.shape[2] + x_indices[:, None, :]
    )
    samples = flat_map[flat_indices]
    weights = y_weights[:, :, None] * x_weights[:, None, :]
    return np.einsum("nxyc,nxy->nc", samples, weights, optimize=False)


def _cubic(distance: np.ndarray) -> np.ndarray:
    distance = np.abs(distance)
    return np.where(
        distance <= 1,
        1.25 * distance**3 - 2.25 * distance**2 + 1,
        np.where(distance < 2, -0.75 * distance**3 + 3.75 * distance**2 - 6 * distance + 3, 0),
    )


