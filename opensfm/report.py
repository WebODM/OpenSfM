import logging
import os
import subprocess
import tempfile

import PIL
from fpdf import FPDF
from opensfm import io
from opensfm.dataset import DataSet
from typing import Any, Dict

logger: logging.Logger = logging.getLogger(__name__)

UNIT_SYSTEMS = {
    "m": {
        "factor": 1.0, 
        "label": "m", 
        "area_label": "km²", 
        "area_factor": 1 / 1e6,
        "gsd_label": "cm", 
        "gsd_factor": 1.0
    },
    "ft": {
        "factor": 1.0 / 0.3048, 
        "label": "ft", 
        "area_label": "ac", 
        "area_factor": (1 / (0.3048 * 0.3048)) / 43560, 
        "gsd_label": "in", 
        "gsd_factor": 1.0 / 2.54
    },
    "US survey foot": {
        "factor": 3937 / 1200, 
        "label": "ft US", 
        "area_label": "ac US", 
        "area_factor": (3937 / 1200)**2 / 43560, 
        "gsd_label": "in", 
        "gsd_factor": 1.0 / 2.54
    },
}


class Report:
    def __init__(self, data: DataSet, stats = None) -> None:
        self.output_path = os.path.join(data.data_path, "stats")
        self.dataset_name = os.path.basename(data.data_path)
        self.io_handler = data.io_handler
        units = data.config.get("report_unit_system", "m")

        if units not in UNIT_SYSTEMS:
            units = "m"

        us = UNIT_SYSTEMS[units]
        self.unit_factor = us["factor"]
        self.unit_label = us["label"]
        self.area_label = us["area_label"]
        self.area_factor = us["area_factor"]
        self.gsd_label = us["gsd_label"]
        self.gsd_factor = us["gsd_factor"]

        self.mapi_light_light_green = [255, 255, 255]
        self.mapi_light_green = [0, 0, 0]
        self.mapi_light_grey = [243, 244, 247]
        self.mapi_dark_grey = [0, 0, 0]
        self.border_color = [100, 100, 100]

        self.pdf = FPDF("P", "mm", "A4")
        self.pdf.add_page()

        self.title_size = 20
        self.h1 = 16
        self.h2 = 13
        self.h3 = 10
        self.text = 10
        self.small_text = 8
        self.margin = 10
        self.cell_height = 7
        self.total_size = 190

        if stats is not None:
            self.stats = stats
        else:
            self.stats = self._read_stats_file("stats.json")
        
        self.version = data.config["report_version"]
        self.name = data.config["report_name"]

    def save_report(self, filename: str) -> None:
        # pyre-fixme[28]: Unexpected keyword argument `dest`.
        bytestring = self.pdf.output(dest="S")
        if isinstance(bytestring, str):
            bytestring = bytestring.encode("utf8")

        with self.io_handler.open_wb(
            os.path.join(self.output_path, filename)
        ) as fwb:
            fwb.write(bytestring)

    def _make_table(self, columns_names, rows, row_header=False, columns_sizes=None) -> None:
        if len(rows) == 0:
            logger.warning("Cannot make table (rows missing)")
            return

        self.pdf.set_font("Helvetica", "B", self.h3)
        self.pdf.set_line_width(0.3)

        if columns_sizes is None:
            columns_sizes = [int(self.total_size / len(rows[0]))] * len(rows[0])
        else:
            columns_sizes = [int(self.total_size) * p for p in columns_sizes]

        if columns_names:
            self.pdf.set_draw_color(*self.border_color)
            self.pdf.set_fill_color(*self.mapi_light_grey)
            for col, size in zip(columns_names, columns_sizes):
                self.pdf.rect(
                    self.pdf.get_x(),
                    self.pdf.get_y(),
                    size,
                    self.cell_height,
                    style="FD",
                )
                self.pdf.set_text_color(*self.mapi_dark_grey)
                self.pdf.cell(size, self.cell_height, col, align="L")
            self.pdf.set_xy(self.margin, self.pdf.get_y() + self.cell_height)

        self.pdf.set_draw_color(*self.border_color)
        self.pdf.set_fill_color(*self.mapi_light_light_green)

        self.pdf.set_font("Helvetica", "", self.h3)

        for row in rows:
            for i, (col, size) in enumerate(zip(row, columns_sizes)):
                if i == 0 and row_header:
                    self.pdf.set_draw_color(*self.border_color)
                    self.pdf.set_fill_color(*self.mapi_light_grey)
                self.pdf.rect(
                    self.pdf.get_x(),
                    self.pdf.get_y(),
                    size,
                    self.cell_height,
                    style="FD",
                )
                self.pdf.set_text_color(*self.mapi_dark_grey)
                if i == 0 and row_header:
                    self.pdf.set_font("Helvetica", "B", self.h3)
                    self.pdf.set_draw_color(*self.border_color)
                    self.pdf.set_fill_color(*self.mapi_light_light_green)
                else:
                    self.pdf.set_font("Helvetica", "", self.h3)
                self.pdf.cell(size, self.cell_height, col, align="L")
            self.pdf.set_xy(self.margin, self.pdf.get_y() + self.cell_height)

    def _read_stats_file(self, filename) -> Dict[str, Any]:
        file_path = os.path.join(self.output_path, filename)
        with self.io_handler.open_rt(file_path) as fin:
            return io.json_load(fin)

    def _get_gcp_stats(self):
        return self.stats.get("gcp_errors", {}).get("details", [])

    def _to_unit(self, value_m: float) -> float:
        """Convert a value in meters to the configured unit."""
        return value_m * self.unit_factor

    def _to_area(self, value_sqm: float) -> float:
        """Convert a value in square meters to the configured area unit."""
        return value_sqm * self.area_factor

    def _to_gsd(self, value_cm: float) -> float:
        """Convert a GSD value in cm to the configured small unit."""
        return value_cm * self.gsd_factor

    def _make_section(self, title: str) -> None:
        self.pdf.set_font("Helvetica", "B", self.h1)
        self.pdf.set_text_color(*self.mapi_dark_grey)
        self.pdf.cell(0, self.margin, title, align="L")
        self.pdf.set_xy(self.margin, self.pdf.get_y() + 1.5 * self.margin)

    def _make_subsection(self, title: str) -> None:
        self.pdf.set_xy(self.margin, self.pdf.get_y() - 0.5 * self.margin)
        self.pdf.set_font("Helvetica", "B", self.h2)
        self.pdf.set_text_color(*self.mapi_dark_grey)
        self.pdf.cell(0, self.margin, title, align="L")
        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin)

    def _make_centered_image(self, image_path: str, desired_height: float) -> None:

        with tempfile.TemporaryDirectory() as tmp_local_dir:
            local_image_path = os.path.join(tmp_local_dir, os.path.basename(image_path))
            with self.io_handler.open_wb(local_image_path) as fwb:
                with self.io_handler.open_rb(image_path) as f:
                    fwb.write(f.read())

            width, height = PIL.Image.open(local_image_path).size
            resized_width = width * desired_height / height
            if resized_width > self.total_size:
                resized_width = self.total_size
                desired_height = height * resized_width / width

            self.pdf.image(
                local_image_path,
                self.pdf.get_x() + self.total_size / 2 - resized_width / 2,
                self.pdf.get_y(),
                h=desired_height,
            )
            self.pdf.set_xy(
                self.margin, self.pdf.get_y() + desired_height + self.margin
            )

    def make_title(self) -> None:
        # title
        self.pdf.set_font("Helvetica", "B", self.title_size)
        self.pdf.set_text_color(*self.mapi_light_green)
        self.pdf.cell(0, self.margin, "Quality Report", align="C")
        self.pdf.set_xy(self.margin, self.title_size)

        # version number
        version = f"version {self.version}" if self.version != "" else ""

        self.pdf.set_font("Helvetica", "", self.small_text)
        self.pdf.set_text_color(*self.mapi_dark_grey)
        self.pdf.cell(
            0, self.margin, f"Processed with {self.name} {version}", align="R"
        )
        self.pdf.set_xy(self.margin, self.pdf.get_y() + 2 * self.margin)

    def make_dataset_summary(self) -> None:
        self._make_section("Dataset Summary")

        rows = [
            ["Date", self.stats["processing_statistics"]["date"]],
            [
                "Area Covered",
                f"{self._to_area(self.stats['processing_statistics']['area']):.6f} {self.area_label}",
            ],
            [
                "Processing Time",
                #f"{self.stats['processing_statistics']['steps_times']['Total Time']:.2f} seconds",
                self.stats['odm_processing_statistics']['total_time_human'],
            ],
            ["Capture Start", self.stats["processing_statistics"]["start_date"]],
            ["Capture End", self.stats["processing_statistics"]["end_date"]],
        ]

        if "crs" in self.stats and "crs_name" in self.stats:
            crs = self.stats["crs"]
            if crs == "":
                crs = "Local"
            elif self.stats["crs_name"] and self.stats["crs_name"] != "unknown":
                crs = f"{self.stats['crs_name']}"
            rows.append(["Coordinate Reference System", crs])

        self._make_table(None, rows, True, columns_sizes=[0.35, 0.65])
        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

    def _has_meaningful_gcp(self) -> bool:
        return (
            self.stats["reconstruction_statistics"]["has_gcp"]
            and "average_error" in self.stats["gcp_errors"]
        )

    def make_processing_summary(self) -> None:
        self._make_section("Processing Summary")

        rec_shots, init_shots = (
            self.stats["reconstruction_statistics"]["reconstructed_shots_count"],
            self.stats["reconstruction_statistics"]["initial_shots_count"],
        )
        rec_points, init_points = (
            self.stats["reconstruction_statistics"]["reconstructed_points_count"],
            self.stats["reconstruction_statistics"]["initial_points_count"],
        )

        geo_string = []
        if self.stats["reconstruction_statistics"]["has_gps"]:
            geo_string.append("GPS")
        if self._has_meaningful_gcp():
            geo_string.append("GCP")
        
        if "align" in self.stats:
            geo_string = ["Alignment"]

        if len(geo_string) == 0:
            geo_string = ["None"]

        ratio_shots = rec_shots / init_shots * 100 if init_shots > 0 else -1
        rows = [
            [
                "Reconstructed Images",
                f"{rec_shots} over {init_shots} shots ({ratio_shots:.1f}%)",
            ],
            [
                "Reconstructed Points (Sparse)",
                f"{rec_points} over {init_points} points ({rec_points/init_points*100:.1f}%)",
            ],
            # [
            #     "Reconstructed Components",
            #     f"{self.stats['reconstruction_statistics']['components']} component",
            # ],
            [
                "Detected Features",
                f"{self.stats['features_statistics']['detected_features']['median']:,} features",
            ],
            [
                "Reconstructed Features",
                f"{self.stats['features_statistics']['reconstructed_features']['median']:,} features",
            ],
            ["Geographic Reference", " and ".join(geo_string)],
        ]

        # Dense (if available)
        if self.stats.get('point_cloud_statistics'):
            if self.stats['point_cloud_statistics'].get('dense'):
                rows.insert(2, [
                    "Reconstructed Points (Dense)",
                    f"{self.stats['point_cloud_statistics']['stats']['statistic'][0]['count']:,} points"
                ])

        # GSD (if available)
        if self.stats['odm_processing_statistics'].get('average_gsd'):
            rows.insert(3, [
                "Average Ground Sampling Distance",
                f"{self._to_gsd(self.stats['odm_processing_statistics']['average_gsd']):.2f} {self.gsd_label}"
            ])
        
        row_gps_gcp = [" / ".join(geo_string) + " errors"]
        geo_errors = []
        
        if not "align" in self.stats:
            if self.stats["reconstruction_statistics"]["has_gps"]:
                geo_errors.append(f"{self._to_unit(self.stats['gps_errors']['average_error']):.2f}")
            if self._has_meaningful_gcp():
                geo_errors.append(f"{self._to_unit(self.stats['gcp_errors']['average_error']):.2f}")
        else:
            geo_errors.append(f"{self._to_unit(self.stats['align']['coarse']['rmse_3d'] + self.stats['align']['fine']['rmse_3d']):.2f}")
        
        if len(geo_errors) > 0:
            row_gps_gcp.append(" / ".join(geo_errors) + f" {self.unit_label}")
            rows.append(row_gps_gcp)

        self._make_table(None, rows, True, columns_sizes=[0.35, 0.65])
        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

        topview_height = 110
        topview_grids = [
            f for f in self.io_handler.ls(self.output_path) if f.startswith("topview")
        ]
        self._make_centered_image(
            os.path.join(self.output_path, topview_grids[0]), topview_height
        )

        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin)

    def make_processing_time_details(self) -> None:
        self._make_section("Processing Time Details")

        columns_names = list(self.stats["processing_statistics"]["steps_times"].keys())
        formatted_floats = []
        for v in self.stats["processing_statistics"]["steps_times"].values():
            formatted_floats.append(f"{v:.2f} sec.")
        rows = [formatted_floats]
        self._make_table(columns_names, rows)
        self.pdf.set_xy(self.margin, self.pdf.get_y() + 2 * self.margin)

    def make_gcp_error_details(self):
        self._make_section("Ground Control Points")
        gcp_stats = self._get_gcp_stats()

        gcp_rows = []
        chk_rows = []

        column_names = ["ID", f"Error X ({self.unit_label})", f"Error Y ({self.unit_label})", f"Error Z ({self.unit_label})"]

        for gcp in gcp_stats:
            if gcp['error'] is None:
                continue
            gcp_id = ''.join(c if ord(c) < 256 else '?' for c in gcp["id"]) # latin-1 only due to fpdf2 requirements
            row = [gcp_id]
            row.append(f"{self._to_unit(gcp['error']['x']):.3f}")
            row.append(f"{self._to_unit(gcp['error']['y']):.3f}")
            row.append(f"{self._to_unit(gcp['error']['z']):.3f}")
            
            if gcp_id.startswith("CHK-"):
                chk_rows.append(row)
            else:
                gcp_rows.append(row)
        
        self._make_table(column_names, gcp_rows)
        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

        if len(chk_rows) > 0:
            self._make_section("Checkpoints")
            self._make_table(column_names, chk_rows)
            self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

    def make_gps_details(self) -> None:
        self._make_section("GPS/GCP/3D Details")

        # GPS
        table_count = 0
        for error_type in ["gps", "gcp", "3d"]:
            stats = self.stats[error_type + "_errors"]
            if error_type == "gcp" and "gcp_only" in stats:
                stats = stats["gcp_only"]
            
            rows = []
            columns_names = [error_type.upper(), "Mean", "Standard Deviation", "RMS Error"]
            if "average_error" not in stats:
                continue
            for comp in ["x", "y", "z"]:
                row = [comp.upper() + f" Error ({self.unit_label})"]
                row.append(f"{self._to_unit(stats['mean'][comp]):.3f}")
                row.append(f"{self._to_unit(stats['std'][comp]):.3f}")
                row.append(f"{self._to_unit(stats['error'][comp]):.3f}")
                rows.append(row)

            rows.append(
                [
                    "Total",
                    "",
                    "",
                    f"{self._to_unit(stats['average_error']):.3f}",
                ]
            )
            self._make_table(columns_names, rows, True)
            self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)
            table_count += 1

        if table_count > 0:
            abs_error_type = "gps" if table_count == 2 else "gcp"

            a_ce90 = self.stats[abs_error_type + "_errors"].get("ce90", 0)
            a_le90 = self.stats[abs_error_type + "_errors"].get("le90", 0)
            r_ce90 = self.stats["3d_errors"].get("ce90", 0)
            r_le90 = self.stats["3d_errors"].get("le90", 0)

            rows = []
            if a_ce90 > 0 and a_le90 > 0:
                rows += [[
                    f"Horizontal CE90 ({self.unit_label})",
                    f"{self._to_unit(a_ce90):.3f}",
                    f"{self._to_unit(r_ce90):.3f}" if r_ce90 > 0 else "-",
                ],[
                    f"Vertical LE90 ({self.unit_label})",
                    f"{self._to_unit(a_le90):.3f}",
                    f"{self._to_unit(r_le90):.3f}" if r_le90 > 0 else "-",
                ]]
            
            if rows:
                if table_count > 2:
                    self.add_page_break()
                self._make_table(["Accuracy", "Absolute", "Relative"], rows, True)
                self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

        # rows = []
        # columns_names = [
        #     "GPS Bias",
        #     "Scale",
        #     "Translation",
        #     "Rotation",
        # ]
        # for camera, params in self.stats["camera_errors"].items():
        #     bias = params["bias"]
        #     s, t, R = bias["scale"], bias["translation"], bias["rotation"]
        #     rows.append(
        #         [
        #             camera,
        #             f"{s:.2f}",
        #             f"{t[0]:.2f}      {t[1]:.2f}      {t[2]:.2f}",
        #             f"{R[0]:.2f}      {R[1]:.2f}      {R[2]:.2f}",
        #         ]
        #     )
        # self._make_table(columns_names, rows)

        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

    def make_align_details(self) -> None:
        self._make_section("Alignment Errors Details")

        # Alignment
        rows = []
        columns_names = ["", "DSM (Coarse)", "ICP (Fine)", "RMS Error"]
        for comp in ["x", "y", "z", "3d"]:
            row = [comp.upper() + f" Error ({self.unit_label})"]
            row.append(f"{self._to_unit(self.stats['align']['coarse']['rmse_' + comp]):.3f}")
            row.append(f"{self._to_unit(self.stats['align']['fine']['rmse_' + comp]):.3f}")
            row.append(f"{self._to_unit(self.stats['align']['coarse']['rmse_' + comp] + self.stats['align']['fine']['rmse_' + comp]):.3f}")
            rows.append(row)

        self._make_table(columns_names, rows)
        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

        dsm_feature_matches = os.path.join(self.output_path, "codem", "dsm_feature_matches.png")
        if os.path.isfile(dsm_feature_matches):
            self._make_centered_image(dsm_feature_matches, 80)

    def make_features_details(self) -> None:
        self._make_section("Features Details")

        heatmap_height = 60
        heatmaps = [
            f for f in self.io_handler.ls(self.output_path) if f.startswith("heatmap") and f.endswith(".png")
        ]
        self._make_centered_image(
            os.path.join(self.output_path, heatmaps[0]), heatmap_height
        )
        if len(heatmaps) > 1:
            logger.warning("Please implement multi-model display")

        columns_names = ["", "Min.", "Max.", "Mean", "Median"]
        rows = []
        for comp in ["detected_features", "reconstructed_features"]:
            row = [comp.replace("_", " ").replace("features", "").capitalize()]
            for t in columns_names[1:]:
                row.append(
                    f"{self.stats['features_statistics'][comp][t.replace('.', '').lower()]:.0f}"
                )
            rows.append(row)
        self._make_table(columns_names, rows, True)

        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin)

    def make_reconstruction_details(self) -> None:
        self._make_section("Reconstruction Details")

        rows = [
            [
                "Average Reprojection Error (norm / pixels / angular)",
                (
                    f"{self.stats['reconstruction_statistics']['reprojection_error_normalized']:.2f} / "
                    f"{self.stats['reconstruction_statistics']['reprojection_error_pixels']:.2f} / "
                    f"{self.stats['reconstruction_statistics']['reprojection_error_angular']:.5f}"
                ),
            ],
            [
                "Average Track Length",
                f"{self.stats['reconstruction_statistics']['average_track_length']:.2f} images",
            ],
            [
                "Average Track Length (> 2)",
                f"{self.stats['reconstruction_statistics']['average_track_length_over_two']:.2f} images",
            ],
        ]
        self._make_table(None, rows, True)
        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 1.5)

        residual_histogram_height = 60
        residual_histogram = [
            f
            for f in self.io_handler.ls(self.output_path)
            if f.startswith("residual_histogram")
        ]
        self._make_centered_image(
            os.path.join(self.output_path, residual_histogram[0]),
            residual_histogram_height,
        )
        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin)

    def make_camera_models_details(self) -> None:
        self._make_section("Camera Models Details")

        for camera, params in self.stats["camera_errors"].items():
            residual_grids = [
                f
                for f in self.io_handler.ls(self.output_path)
                if f.startswith("residuals_" + str(camera.replace("/", "_"))) and f.endswith(".png")
            ]
            if not residual_grids:
                continue

            initial = params["initial_values"]
            optimized = params["optimized_values"]
            names = [""] + list(initial.keys())
            names = [n if n != "aspect_ratio" else "aspect" for n in names]

            rows = []
            rows.append(["Initial"] + [f"{x:.4f}" for x in initial.values()])
            rows.append(["Optimized"] + [f"{x:.4f}" for x in optimized.values()])

            self._make_subsection(camera)
            self._make_table(names, rows, True)
            self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

            residual_grid_height = 100
            self._make_centered_image(
                os.path.join(self.output_path, residual_grids[0]), residual_grid_height
            )

    def make_rig_cameras_details(self) -> None:
        if len(self.stats["rig_errors"]) == 0:
            return

        self._make_section("Rig Cameras Details")

        columns_names = [
            "Translation X",
            "Translation Y",
            "Translation Z",
            "Rotation X",
            "Rotation Y",
            "Rotation Z",
        ]
        for rig_camera_id, params in self.stats["rig_errors"].items():
            initial = params["initial_values"]
            optimized = params["optimized_values"]

            rows = []
            r_init, t_init = initial["rotation"], initial["translation"]
            r_opt, t_opt = optimized["rotation"], optimized["translation"]
            rows.append(
                [
                    f"{self._to_unit(t_init[0]):.4f} {self.unit_label}",
                    f"{self._to_unit(t_init[1]):.4f} {self.unit_label}",
                    f"{self._to_unit(t_init[2]):.4f} {self.unit_label}",
                    f"{r_init[0]:.4f}",
                    f"{r_init[1]:.4f}",
                    f"{r_init[2]:.4f}",
                ]
            )
            rows.append(
                [
                    f"{self._to_unit(t_opt[0]):.4f} {self.unit_label}",
                    f"{self._to_unit(t_opt[1]):.4f} {self.unit_label}",
                    f"{self._to_unit(t_opt[2]):.4f} {self.unit_label}",
                    f"{r_opt[0]:.4f}",
                    f"{r_opt[1]:.4f}",
                    f"{r_opt[2]:.4f}",
                ]
            )

            self._make_subsection(rig_camera_id)
            self._make_table(columns_names, rows)
            self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)

    def make_tracks_details(self) -> None:
        self._make_section("Tracks Details")
        matchgraph_height = 80
        matchgraph = [
            f
            for f in self.io_handler.ls(self.output_path)
            if f.startswith("matchgraph") and f.endswith(".png")
        ]
        self._make_centered_image(
            os.path.join(self.output_path, matchgraph[0]), matchgraph_height
        )

        histogram = self.stats["reconstruction_statistics"]["histogram_track_length"]
        start_length, end_length = 2, 10
        row_length = ["Length"]
        for length, _ in sorted(histogram.items(), key=lambda x: int(x[0])):
            if int(length) < start_length or int(length) > end_length:
                continue
            row_length.append(length)
        row_count = ["Count"]
        for length, count in sorted(histogram.items(), key=lambda x: int(x[0])):
            if int(length) < start_length or int(length) > end_length:
                continue
            row_count.append(f"{count}")

        self._make_table(None, [row_length, row_count], True)

        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin)

    def add_page_break(self) -> None:
        self.pdf.add_page("P")

    def make_survey_data(self):
        self._make_section("Survey Data")

        self._make_centered_image(
            os.path.join(self.output_path, "overlap.png"), 90
        )
        self._make_centered_image(
            os.path.join(self.output_path, "overlap_diagram_legend.png"), 3
        )

        self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin / 2)


    def _add_image_label(self, text):
        self.pdf.set_font_size(self.small_text)
        self.pdf.text(self.pdf.get_x() + self.total_size / 2 - self.pdf.get_string_width(text) / 2, self.pdf.get_y() - 5, text)


    def make_preview(self):
        ortho = os.path.join(self.output_path, "ortho.png")
        dsm = os.path.join(self.output_path, "dsm.png")
        dtm = os.path.join(self.output_path, "dtm.png")
        count = 0

        if os.path.isfile(ortho) or os.path.isfile(dsm):
            self._make_section("Previews")
            
            if os.path.isfile(ortho):
                self._make_centered_image(
                    os.path.join(self.output_path, ortho), 110
                )
                self._add_image_label("Orthophoto")
                count += 1

            if os.path.isfile(dsm) and self.stats.get('dsm_statistics'):
                self._make_centered_image(
                    os.path.join(self.output_path, dsm), 110
                )
                self._add_image_label("Digital Surface Model")

                self._make_centered_image(
                    os.path.join(self.output_path, "dsm_gradient.png"), 4
                )
                self.pdf.set_font_size(self.small_text)
                min_text = "{:,.2f} {}".format(self._to_unit(self.stats['dsm_statistics']['min']), self.unit_label)
                max_text = "{:,.2f} {}".format(self._to_unit(self.stats['dsm_statistics']['max']), self.unit_label)
                self.pdf.text(self.pdf.get_x() + 40, self.pdf.get_y() - 5, min_text)
                self.pdf.text(self.pdf.get_x() + 40 + 110.5 - self.pdf.get_string_width(max_text), self.pdf.get_y() - 5, max_text)
                count += 1

            if os.path.isfile(dtm) and self.stats.get('dtm_statistics'):
                if count >= 2:
                    self.add_page_break()

                self._make_centered_image(
                    os.path.join(self.output_path, dtm), 110
                )
                self._add_image_label("Digital Terrain Model")

                self._make_centered_image(
                    os.path.join(self.output_path, "dsm_gradient.png"), 4
                )
                self.pdf.set_font_size(self.small_text)
                min_text = "{:,.2f} {}".format(self._to_unit(self.stats['dtm_statistics']['min']), self.unit_label)
                max_text = "{:,.2f} {}".format(self._to_unit(self.stats['dtm_statistics']['max']), self.unit_label)
                self.pdf.text(self.pdf.get_x() + 40, self.pdf.get_y() - 5, min_text)
                self.pdf.text(self.pdf.get_x() + 40 + 110.5 - self.pdf.get_string_width(max_text), self.pdf.get_y() - 5, max_text)

            self.pdf.set_xy(self.margin, self.pdf.get_y() + self.margin)

            return True

    def generate_report(self) -> None:
        self.make_title()
        self.make_dataset_summary()
        self.make_processing_summary()
        self.add_page_break()

        if self.make_preview():
            self.add_page_break()

        if os.path.isfile(os.path.join(self.output_path, "overlap.png")):
            self.make_survey_data()
        
        if "align" not in self.stats:
            self.make_gps_details()

            gcp_stats = self._get_gcp_stats()
            if len(gcp_stats) > 0:
                self.make_gcp_error_details()
        else:
            self.make_align_details()

        self.add_page_break()

        self.make_features_details()
        self.make_reconstruction_details()
        self.add_page_break()

        self.make_tracks_details()
        self.make_camera_models_details()
        #self.make_rig_cameras_details()
