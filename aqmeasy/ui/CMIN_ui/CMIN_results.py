"""Results display for CMIN execution"""
import os

from PySide6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QGroupBox,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QTabWidget,
    QComboBox,
    QToolButton,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QPainter, QColor, QPen

try:
    from PySide6.QtCharts import (
        QBarCategoryAxis,
        QBarSeries,
        QBarSet,
        QChart,
        QChartView,
        QLineSeries,
        QValueAxis,
    )
    HAS_QT_CHARTS = True
except Exception:
    HAS_QT_CHARTS = False
    QBarCategoryAxis = None
    QBarSeries = None
    QBarSet = None
    QChart = None
    QChartView = None
    QLineSeries = None
    QValueAxis = None


class ResultsPanel(QWidget):
    """Widget to display CMIN execution results"""
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        results_group = QGroupBox("Results Summary")
        results_layout = QVBoxLayout()
        results_group.setLayout(results_layout)
        layout.addWidget(results_group)

        stats_layout = QHBoxLayout()
        results_layout.addLayout(stats_layout)

        input_box = QVBoxLayout()
        self.input_count_label = QLabel("0")
        self.input_count_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.input_count_label.setStyleSheet("font-size: 24px; font-weight: bold; color: #3498db;")
        input_label = QLabel("Input Files")
        input_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        input_box.addWidget(self.input_count_label)
        input_box.addWidget(input_label)
        stats_layout.addLayout(input_box)

        output_box = QVBoxLayout()
        self.output_count_label = QLabel("0")
        self.output_count_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.output_count_label.setStyleSheet("font-size: 24px; font-weight: bold; color: #2ecc71;")
        output_label = QLabel("Output Files")
        output_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        output_box.addWidget(self.output_count_label)
        output_box.addWidget(output_label)
        stats_layout.addLayout(output_box)

        file_eliminated_box = QVBoxLayout()
        self.file_eliminated_label = QLabel("0")
        self.file_eliminated_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.file_eliminated_label.setStyleSheet("font-size: 24px; font-weight: bold; color: #e67e22;")
        file_eliminated_text = QLabel("Filtered to Zero")
        file_eliminated_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        file_eliminated_box.addWidget(self.file_eliminated_label)
        file_eliminated_box.addWidget(file_eliminated_text)
        stats_layout.addLayout(file_eliminated_box)

        failed_files_box = QVBoxLayout()
        self.failed_files_label = QLabel("0")
        self.failed_files_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.failed_files_label.setStyleSheet("font-size: 24px; font-weight: bold; color: #c0392b;")
        failed_files_text = QLabel("Execution Failed")
        failed_files_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        failed_files_box.addWidget(self.failed_files_label)
        failed_files_box.addWidget(failed_files_text)
        stats_layout.addLayout(failed_files_box)

        conformer_eliminated_box = QVBoxLayout()
        self.conformer_eliminated_label = QLabel("0")
        self.conformer_eliminated_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.conformer_eliminated_label.setStyleSheet("font-size: 24px; font-weight: bold; color: #e74c3c;")
        conformer_eliminated_text = QLabel("Conformer-level Eliminated")
        conformer_eliminated_text.setAlignment(Qt.AlignmentFlag.AlignCenter)
        conformer_eliminated_box.addWidget(self.conformer_eliminated_label)
        conformer_eliminated_box.addWidget(conformer_eliminated_text)
        stats_layout.addLayout(conformer_eliminated_box)

        self.output_dir_label = QLabel("Output directory: Not set")
        self.output_dir_label.setWordWrap(True)
        results_layout.addWidget(self.output_dir_label)

        self.conformer_totals_label = QLabel("Conformers: input 0 | output 0")
        self.conformer_totals_label.setWordWrap(True)
        results_layout.addWidget(self.conformer_totals_label)

        self.elimination_definition_label = QLabel(
            "Filtered to zero counts only files that ran successfully but kept zero conformers."
        )
        self.elimination_definition_label.setWordWrap(True)
        self.elimination_definition_label.setStyleSheet("color: #566573;")
        results_layout.addWidget(self.elimination_definition_label)

        self.warning_label = QLabel("")
        self.warning_label.setWordWrap(True)
        self.warning_label.setStyleSheet("color: #b9770e;")
        results_layout.addWidget(self.warning_label)

        self.details_toggle = QToolButton()
        self.details_toggle.setText("Show detailed results")
        self.details_toggle.setCheckable(True)
        self.details_toggle.setChecked(False)
        self.details_toggle.toggled.connect(self._toggle_details)
        layout.addWidget(self.details_toggle)

        self.results_tabs = QTabWidget()
        layout.addWidget(self.results_tabs)
        self.results_tabs.setVisible(False)

        self.results_tabs.setSizePolicy(
            self.results_tabs.sizePolicy().horizontalPolicy(),
            self.results_tabs.sizePolicy().verticalPolicy(),
        )

        table_tab = QWidget()
        table_layout = QVBoxLayout()
        table_tab.setLayout(table_layout)

        selector_row = QHBoxLayout()
        selector_row.addWidget(QLabel("Selected file:"))
        self.file_selector = QComboBox()
        self.file_selector.currentIndexChanged.connect(self._on_file_selection_changed)
        selector_row.addWidget(self.file_selector)
        table_layout.addLayout(selector_row)

        self.file_table = QTableWidget(0, 10)
        self.file_table.setHorizontalHeaderLabels([
            "Input File",
            "Execution",
            "Filtering",
            "Input Conf",
            "Output Conf",
            "Eliminated Conf",
            "Output File",
            "Energy Min",
            "Energy Max",
            "Structures",
        ])
        self.file_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.file_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.file_table.horizontalHeader().setSectionResizeMode(6, QHeaderView.ResizeMode.Stretch)
        self.file_table.setAlternatingRowColors(True)
        table_layout.addWidget(self.file_table)

        self.structure_table = QTableWidget(0, 4)
        self.structure_table.setHorizontalHeaderLabels([
            "Rank",
            "Structure",
            "Energy",
            "RMSD to Lowest (A)",
        ])
        self.structure_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.structure_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.structure_table.setAlternatingRowColors(True)
        table_layout.addWidget(self.structure_table)

        self.results_tabs.addTab(table_tab, "Per-file Results")

        energy_tab = QWidget()
        energy_outer_layout = QVBoxLayout()
        energy_layout = QHBoxLayout()

        graph_scope_row = QHBoxLayout()
        graph_scope_row.addWidget(QLabel("Graph scope:"))
        self.graph_scope_selector = QComboBox()
        self.graph_scope_selector.addItems(["Selected file", "All files"])
        self.graph_scope_selector.currentIndexChanged.connect(self._refresh_charts)
        graph_scope_row.addWidget(self.graph_scope_selector)
        graph_scope_row.addStretch(1)
        energy_outer_layout.addLayout(graph_scope_row)
        energy_outer_layout.addLayout(energy_layout)
        energy_tab.setLayout(energy_outer_layout)

        charts_enabled = HAS_QT_CHARTS and all(
            cls is not None
            for cls in (QChart, QChartView, QBarCategoryAxis, QBarSeries, QBarSet, QLineSeries, QValueAxis)
        )

        if charts_enabled:
            assert QChart is not None
            assert QChartView is not None
            self.hist_chart = QChart()
            self.hist_chart.setTitle("Energy Distribution")
            self.hist_chart.legend().setVisible(False)
            self.hist_chart_view = QChartView(self.hist_chart)
            self.hist_chart_view.setRenderHint(QPainter.RenderHint.Antialiasing)
            self.hist_chart_view.setMinimumHeight(300)
            energy_layout.addWidget(self.hist_chart_view)

            self.profile_chart = QChart()
            self.profile_chart.setTitle("Per-file Energy Profiles")
            self.profile_chart.legend().setVisible(True)
            self.profile_chart_view = QChartView(self.profile_chart)
            self.profile_chart_view.setRenderHint(QPainter.RenderHint.Antialiasing)
            self.profile_chart_view.setMinimumHeight(300)
            energy_layout.addWidget(self.profile_chart_view)
            energy_layout.setStretch(0, 1)
            energy_layout.setStretch(1, 1)
        else:
            chart_missing = QLabel("QtCharts module is unavailable in this environment.")
            chart_missing.setAlignment(Qt.AlignmentFlag.AlignCenter)
            energy_layout.addWidget(chart_missing)

        self.results_tabs.addTab(energy_tab, "Energy Graphs")

        self.charts_enabled = charts_enabled
        self._per_file_by_name = {}
        self._latest_energy_values = []
        self._latest_per_file_data = []
        self.hide()

    def update_results(self, results):
        """Update display with results dictionary"""
        self.input_count_label.setText(str(results.get('input_count', 0)))
        self.output_count_label.setText(str(results.get('output_count', 0)))
        self.file_eliminated_label.setText(str(results.get('file_level_eliminated', results.get('eliminated_count', 0))))
        self.failed_files_label.setText(str(results.get('failed_file_count', 0)))
        self.conformer_eliminated_label.setText(str(results.get('conformer_level_eliminated', 0)))

        if results.get('output_dir'):
            self.output_dir_label.setText(f"Output directory: {results['output_dir']}/CMIN")

        input_conformers = results.get('input_conformer_count', 0)
        output_conformers = results.get('output_conformer_count', 0)
        self.conformer_totals_label.setText(
            f"Conformers: input {input_conformers} | output {output_conformers}"
        )

        warnings = results.get('warnings', [])
        if warnings:
            self.warning_label.setText("Warnings: " + " | ".join(warnings))
        else:
            self.warning_label.setText("")

        self._populate_file_table(results.get('per_file_data', []))
        self._latest_energy_values = list(results.get('energy_values', []))
        self._latest_per_file_data = list(results.get('per_file_data', []))
        self._refresh_charts()

        if results.get('per_file_data'):
            self.file_selector.setCurrentIndex(0)
            self._on_file_selection_changed(0)
            self.details_toggle.setChecked(True)

        self.show()

    def _populate_file_table(self, per_file_data):
        """Populate detailed per-file results table."""
        self.file_table.setRowCount(len(per_file_data))
        self._per_file_by_name = {}

        self.file_selector.blockSignals(True)
        self.file_selector.clear()

        for row, item in enumerate(per_file_data):
            output_files = item.get('matched_output_files', [])
            output_display = output_files[0] if output_files else "None"
            min_energy = item.get('min_energy')
            max_energy = item.get('max_energy')
            structures = item.get('structures', [])
            input_name = item.get('input_name', '')
            self._per_file_by_name[input_name] = item
            self.file_selector.addItem(input_name)

            cells = [
                input_name,
                self._format_outcome_text(item.get('execution_status', 'unknown')),
                self._format_outcome_text(item.get('filter_outcome', 'unknown')),
                self._format_optional_int(item.get('input_conformers')),
                self._format_optional_int(item.get('output_conformers')),
                self._format_optional_int(item.get('eliminated_conformers')),
                output_display,
                self._format_optional_float(min_energy),
                self._format_optional_float(max_energy),
                str(len(structures)),
            ]

            for column, value in enumerate(cells):
                table_item = QTableWidgetItem(value)
                if column in (3, 4, 5, 9):
                    table_item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.file_table.setItem(row, column, table_item)

        self.file_selector.blockSignals(False)

    def _populate_structure_table(self, structures):
        """Populate per-structure table for selected file."""
        self.structure_table.setRowCount(len(structures))
        for row, structure in enumerate(structures):
            cells = [
                self._format_optional_int(structure.get('rank')),
                structure.get('name', f"Conformer {row + 1}"),
                self._format_optional_float(structure.get('energy')),
                self._format_optional_float(structure.get('rmsd_to_best')),
            ]
            for col, value in enumerate(cells):
                item = QTableWidgetItem(value)
                if col in (0, 2, 3):
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                self.structure_table.setItem(row, col, item)

    def _on_file_selection_changed(self, index):
        """Update structure table for selected file record."""
        if index < 0:
            self.structure_table.setRowCount(0)
            self._refresh_charts()
            return
        selected_name = self.file_selector.currentText()
        selected_record = self._per_file_by_name.get(selected_name, {})
        self._populate_structure_table(selected_record.get('structures', []))
        self._refresh_charts()

    def _toggle_details(self, checked):
        """Collapse/expand detailed tabs to keep main window compact."""
        self.results_tabs.setVisible(checked)
        self.details_toggle.setText("Hide detailed results" if checked else "Show detailed results")

    def _refresh_charts(self):
        """Refresh charts using either selected-file or global data scope."""
        if not self.charts_enabled:
            return

        selected_name = self.file_selector.currentText().strip()
        selected_record = self._per_file_by_name.get(selected_name, {})
        selected_structures = selected_record.get('structures', [])
        selected_energies = [
            structure.get('energy')
            for structure in selected_structures
            if structure.get('energy') is not None
        ]

        scope_all = self.graph_scope_selector.currentText() == "All files"
        if scope_all or not selected_energies:
            energy_values = list(self._latest_energy_values)
            per_file_data = list(self._latest_per_file_data)
            chart_title_suffix = " (all files)"
        else:
            energy_values = list(selected_energies)
            per_file_data = [selected_record]
            chart_title_suffix = f" ({selected_name})"

        self._update_charts(energy_values, per_file_data, chart_title_suffix)

    def _update_charts(self, energy_values, per_file_data, title_suffix=""):
        """Render readable summary charts for energies and per-file outcomes."""
        if not self.charts_enabled:
            return

        assert QBarSet is not None
        assert QBarSeries is not None
        assert QBarCategoryAxis is not None
        assert QLineSeries is not None
        assert QValueAxis is not None

        self._reset_chart(self.hist_chart)
        self.hist_chart.legend().setVisible(True)

        if energy_values:
            sorted_energies = sorted(float(value) for value in energy_values)
            reference = sorted_energies[0]
            ladder_series = QLineSeries()
            ladder_series.setName("Delta E from minimum")
            line_pen = QPen(QColor("#1f618d"))
            line_pen.setWidth(3)
            ladder_series.setPen(line_pen)
            for rank, energy in enumerate(sorted_energies, start=1):
                ladder_series.append(float(rank), float(energy - reference))

            self.hist_chart.addSeries(ladder_series)
            axis_x = QValueAxis()
            axis_x.setTitleText("Conformer rank")
            axis_x.setRange(1.0, float(max(1, len(sorted_energies))))
            axis_x.setLabelFormat("%d")

            axis_y = QValueAxis()
            axis_y.setTitleText("Delta E")
            max_delta = max(0.0, max(sorted_energies) - reference)
            axis_y.setRange(0.0, max_delta + 0.1)

            self.hist_chart.addAxis(axis_x, Qt.AlignmentFlag.AlignBottom)
            self.hist_chart.addAxis(axis_y, Qt.AlignmentFlag.AlignLeft)
            ladder_series.attachAxis(axis_x)
            ladder_series.attachAxis(axis_y)
            self.hist_chart.setTitle(f"Relative Energy Ladder{title_suffix}")
        else:
            self.hist_chart.legend().setVisible(False)
            self.hist_chart.setTitle(f"Relative Energy Ladder{title_suffix} (no energies found)")

        self._reset_chart(self.profile_chart)
        self.profile_chart.legend().setVisible(True)

        output_series = QBarSet("Output")
        eliminated_series = QBarSet("Eliminated")
        output_series.setColor(QColor("#17a589"))
        eliminated_series.setColor(QColor("#d35400"))

        categories = []
        max_count = 0
        for index, item in enumerate(per_file_data, start=1):
            input_name = item.get('input_name', f"File {index}")
            basename = os.path.basename(input_name)
            if len(basename) > 16:
                basename = basename[:13] + "..."

            output_count = int(item.get('output_conformers', 0) or 0)
            eliminated_count = int(item.get('eliminated_conformers', 0) or 0)
            max_count = max(max_count, output_count, eliminated_count)

            categories.append(basename)
            output_series.append(output_count)
            eliminated_series.append(eliminated_count)

        if categories:
            bar_series = QBarSeries()
            bar_series.append(output_series)
            bar_series.append(eliminated_series)
            self.profile_chart.addSeries(bar_series)

            axis_x = QBarCategoryAxis()
            axis_x.append(categories)
            axis_x.setLabelsAngle(-20)

            axis_y = QValueAxis()
            axis_y.setTitleText("Conformer count")
            axis_y.setRange(0.0, float(max_count + 1))
            axis_y.setLabelFormat("%d")

            self.profile_chart.addAxis(axis_x, Qt.AlignmentFlag.AlignBottom)
            self.profile_chart.addAxis(axis_y, Qt.AlignmentFlag.AlignLeft)
            bar_series.attachAxis(axis_x)
            bar_series.attachAxis(axis_y)
            self.profile_chart.setTitle(f"Per-file Filtering Outcome{title_suffix}")
        else:
            self.profile_chart.legend().setVisible(False)
            self.profile_chart.setTitle(f"Per-file Filtering Outcome{title_suffix} (no files found)")

    def append_log(self, message):
        """Compatibility no-op; execution log was removed from UI."""
        _ = message

    def _reset_chart(self, chart):
        """Remove prior series and axes so refreshed charts keep one clean axis set."""
        chart.removeAllSeries()
        for axis in list(chart.axes()):
            chart.removeAxis(axis)

    def clear(self):
        """Reset results display"""
        self.input_count_label.setText("0")
        self.output_count_label.setText("0")
        self.file_eliminated_label.setText("0")
        self.failed_files_label.setText("0")
        self.conformer_eliminated_label.setText("0")
        self.output_dir_label.setText("Output directory: Not set")
        self.conformer_totals_label.setText("Conformers: input 0 | output 0")
        self.warning_label.setText("")
        self.file_table.setRowCount(0)
        self.structure_table.setRowCount(0)
        self.file_selector.blockSignals(True)
        self.file_selector.clear()
        self.file_selector.blockSignals(False)
        self._per_file_by_name = {}
        self._latest_energy_values = []
        self._latest_per_file_data = []
        self.details_toggle.setChecked(False)
        self.details_toggle.setText("Show detailed results")

        if self.charts_enabled:
            assert hasattr(self, 'hist_chart')
            assert hasattr(self, 'profile_chart')
            self._reset_chart(self.hist_chart)
            self.hist_chart.setTitle("Relative Energy Ladder")
            self._reset_chart(self.profile_chart)
            self.profile_chart.setTitle("Per-file Filtering Outcome")

        self.hide()

    def _format_optional_int(self, value):
        if value is None:
            return "N/A"
        return str(value)

    def _format_optional_float(self, value):
        if value is None:
            return "N/A"
        return f"{float(value):.8f}"

    def _format_outcome_text(self, value):
        return str(value or "unknown").replace('_', ' ').title()