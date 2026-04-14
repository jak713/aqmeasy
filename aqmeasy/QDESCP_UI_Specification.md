# QDESCP UI Specification
## PySide6 Implementation Guide for a Coding Agent

---

## 1. Architectural Overview

### 1.1 Window Lifecycle

QDESCP is a **standalone, independent window** launched from the main hub. It is not embedded in, docked to, or parented as a widget inside the main window.

**Launch sequence:**
1. User clicks the `QDESCP` button on the main hub window.
2. The hub's slot instantiates `QDESCPWindow` and calls `.show()`.
3. `QDESCPWindow` opens as a fully independent, resizable OS window.
4. The hub remains open and usable. Both windows are alive simultaneously.
5. Closing `QDESCPWindow` destroys the instance; the hub is unaffected.

```python
# In the main hub — the only coupling point
from qdescp_window import QDESCPWindow

class MainHub(QMainWindow):
    def __init__(self):
        super().__init__()
        self._qdescp_window = None  # hold reference to prevent GC

    def open_qdescp(self):
        if self._qdescp_window is None or not self._qdescp_window.isVisible():
            self._qdescp_window = QDESCPWindow()
            self._qdescp_window.show()
        else:
            self._qdescp_window.raise_()
            self._qdescp_window.activateWindow()
```

> **Agent instruction:** The `if not isVisible()` guard prevents duplicate windows if the user clicks the button twice. Do not use `exec_()` — that would block the hub.

### 1.2 Class Hierarchy

```
QDESCPWindow(QMainWindow)
│
├── central_widget: QWidget
│   └── main_layout: QVBoxLayout
│       └── tab_widget: QTabWidget
│           ├── Tab 0: DescriptorsTab(QWidget)
│           └── Tab 1: NMRTab(QWidget)
│
└── status_bar: QStatusBar  ← set via self.setStatusBar(...)
```

Each tab is a **self-contained QWidget subclass** in its own file. `QDESCPWindow` imports and instantiates them — it does not own any domain logic.

### 1.3 File Structure

```
qdescp/
├── qdescp_window.py        ← QDESCPWindow: window shell, tab container
├── descriptors_tab.py      ← DescriptorsTab widget
├── nmr_tab.py              ← NMRTab widget
└── workers.py              ← QThread workers (MCS, xTB run, NMR run)
```

---

## 2. QDESCPWindow — The Window Shell

`QDESCPWindow` is responsible for exactly three things:
1. Setting the window title, minimum size, and icon.
2. Creating the `QTabWidget` and adding the two tab widgets.
3. Owning the `QStatusBar` and exposing a `set_status(msg: str)` method that both tabs can call.

It owns **no parameter widgets, no run buttons, and no domain logic**.

```python
class QDESCPWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("QDESCP")
        self.setMinimumSize(720, 600)

        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)

        self.tab_widget = QTabWidget()
        self.descriptors_tab = DescriptorsTab(status_callback=self.set_status)
        self.nmr_tab = NMRTab(status_callback=self.set_status)
        self.tab_widget.addTab(self.descriptors_tab, "Descriptors")
        self.tab_widget.addTab(self.nmr_tab, "NMR")

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.tab_widget)
        self.setCentralWidget(central)

    def set_status(self, message: str):
        self.status_bar.showMessage(message)
```

> **Agent instruction:** Both tab classes receive `status_callback` as a constructor argument (a `Callable[[str], None]`). They call it to push messages to the shared status bar. They do **not** hold a reference to `QDESCPWindow` itself — this keeps tabs fully decoupled from the shell.

---

## 3. Tab Architecture — Shared Rules

These rules apply to **both** `DescriptorsTab` and `NMRTab`:

- Each tab has its own **Run button**, placed at the bottom of its own layout.
- Each tab has its own **progress bar** (`QProgressBar`), directly above its Run button.
- Each tab manages its own **QThread worker**. It starts the worker on Run, and cancels it if the tab is closed or Run is clicked again while running.
- The Run button text toggles: `"Run Descriptors"` → `"Cancel"` while running, then back.
- All heavy computation (file parsing, xTB calls, NMR fitting) runs in a worker thread. **The main thread must never block.**
- Input validation fires when Run is clicked, before the worker is started. Errors are shown via `QMessageBox.warning(self, "Input Error", "...")`.
- After a successful run, the tab calls `self.status_callback("QDESCP descriptors complete.")` or equivalent.

### Tab Layout Template

Every tab follows this vertical structure:

```
┌─────────────────────────────────────────┐
│  [GroupBox: Input]                      │
│  [GroupBox: Settings]                   │
│  [GroupBox: Advanced]  ← collapsible    │
│  ──────────────────────────────────────  │
│  [QProgressBar]                         │
│  [Run Button]                           │
└─────────────────────────────────────────┘
```

Wrap the GroupBoxes in a `QScrollArea` so the window remains usable if the advanced section is expanded. The progress bar and run button live **outside** the scroll area, pinned to the bottom.

---

## 4. Descriptors Tab (`DescriptorsTab`)

### 4.1 Input GroupBox

**Widget: File List**
- Label: `"Input files"`
- Widget: `QListWidget` with `setSelectionMode(QAbstractItemView.ExtendedSelection)`
- Height: minimum 80px, expandable
- Accepts drag-and-drop: `.sdf`, `.xyz`, `.pdb`, `.csv`
- Two buttons beside the list (vertical `QVBoxLayout` to the right):
  - `"Add Files…"` → `QFileDialog.getOpenFileNames` filtered to `SDF / XYZ / PDB / CSV (*.sdf *.xyz *.pdb *.csv)`
  - `"Remove Selected"` → removes highlighted rows

**Format detection logic (runs immediately when files are added):**
```python
def _detect_mode(self, paths: list[str]) -> str:
    """Returns 'csv' if any path ends in .csv, else 'direct'."""
    return 'csv' if any(p.lower().endswith('.csv') for p in paths) else 'direct'
```

When mode is `'csv'`, display an inline info label beneath the file list:
> ℹ️ CSV input detected — QDESCP will run via the CSEARCH conformer pipeline.

When mode is `'direct'`, display:
> ℹ️ Direct input — structures will be passed to xTB without conformer generation.

Use a `QLabel` with word wrap for this. Update it every time the file list changes.

### 4.2 Settings GroupBox

| Parameter | Widget | Default | Notes |
|---|---|---|---|
| `--qdescp_temp` | `QDoubleSpinBox` | 300.0 | Range 0–9999, suffix ` K` |
| `--qdescp_solvent` | `QComboBox` + manual entry | `""` (gas phase) | Pre-populate with common solvents; allow freetext |
| `--boltz` | `QCheckBox` | ✅ checked | Label: `"Boltzmann weighting"` |
| `--xtb_opt` | `QCheckBox` | ✅ checked | Label: `"Optimise geometry with xTB"` |
| `--qdescp_atoms` | `QLineEdit` | `""` | Placeholder: `e.g. P, C=O, [NH2]` |

**`qdescp_atoms` input note:** Add a `(?)` `QToolButton` beside the field that shows a `QMessageBox` explaining: *"Enter atom symbols (e.g. P) or SMARTS patterns (e.g. C=O) separated by commas. Descriptors will be computed for matching atoms."*

### 4.3 Advanced GroupBox (Collapsible)

Implement as a `QGroupBox` with `setCheckable(True)` and `setChecked(False)` so it collapses by default. Its child widget is hidden when unchecked.

```python
advanced_box = QGroupBox("Advanced Options")
advanced_box.setCheckable(True)
advanced_box.setChecked(False)
# Add parameter widgets as children — they are auto-hidden when unchecked
```

| Parameter | Widget | Default |
|---|---|---|
| `--nprocs` | `QSpinBox` | 1, range 1–64 |
| `--qdescp_csv` | `QLineEdit` | `"QDESCP_parameters"` |
| `--qdescp_cclib` | `QCheckBox` | ☐ unchecked |
| `--qdescp_ext` | `QCheckBox` | ☐ unchecked |

### 4.4 Run Button (Descriptors)

- Text: `"Run Descriptors"`
- Object name: `"run_descriptors_btn"` (for testing)
- Validation before starting worker:
  1. File list is not empty.
  2. If `qdescp_atoms` is non-empty, it is a valid comma-separated list of non-empty strings.
- On validation pass: start `DescriptorsWorker`, toggle button to `"Cancel"`.
- On worker `finished` signal: reset button to `"Run Descriptors"`, set progress to 100%.
- On worker `error` signal: reset button, show `QMessageBox.critical`.

---

## 5. NMR Tab (`NMRTab`)

### 5.1 Input GroupBox

Same file list widget as the Descriptors tab, but filtered to `.sdf` and `.xyz` only:
`"SDF / XYZ (*.sdf *.xyz)"`

No mode detection label needed — NMR always operates on direct structures.

### 5.2 Settings GroupBox

| Parameter | Widget | Default | Notes |
|---|---|---|---|
| `--nmr_nucleus` | `QComboBox` | `"H"` | Options: `H`, `C`, `N`, `P`, `F` |
| `--nmr_exp_data` | `QPushButton` + `QLabel` | `""` | Button: `"Browse…"`, label shows filename |
| `--boltz` | `QCheckBox` | ✅ checked | Label: `"Boltzmann weighting"` |
| `--xtb_opt` | `QCheckBox` | ✅ checked | Label: `"Optimise geometry with xTB"` |

### 5.3 Slope / Intercept Table

This is the critical widget for NMR. Use a `QTableWidget` instead of parallel `QLineEdit` fields to enforce row-parity automatically.

```
┌──────────────────────────────────────────┐
│  Scaling parameters  [+ Add] [− Remove]  │
├──────────────────┬───────────┬───────────┤
│  Atom type       │  Slope    │ Intercept │
├──────────────────┼───────────┼───────────┤
│  H               │  1.0491   │  -0.2248  │
│  C               │  1.0195   │  -0.1065  │
└──────────────────┴───────────┴───────────┘
```

- Column 0: `QComboBox` delegate with atom symbols (`H`, `C`, `N`, `P`, `F`)
- Columns 1 & 2: `QDoubleSpinBox` delegate, range −999 to 999, 4 decimal places
- Pre-populate with one row per nucleus in `nmr_nucleus` selection on tab load
- `"+ Add"` appends a new row; `"− Remove"` deletes the selected row
- On Run, extract as parallel lists: `nmr_atoms`, `nmr_slope`, `nmr_intercept`

> **Agent instruction:** Use `QTableWidget` delegates rather than raw item text so values are always valid floats. The `itemDelegateForColumn` API accepts a `QStyledItemDelegate` subclass — implement `createEditor` to return the appropriate spin box.

### 5.4 Advanced GroupBox (Collapsible)

Same collapsible `QGroupBox` pattern as Descriptors.

| Parameter | Widget | Default |
|---|---|---|
| `--nprocs` | `QSpinBox` | 1, range 1–64 |
| `--nmr_online` | `QCheckBox` | ☐ unchecked |

### 5.5 Run Button (NMR)

- Text: `"Run NMR"`
- Object name: `"run_nmr_btn"`
- Validation:
  1. File list is not empty.
  2. Slope/intercept table has at least one row.
  3. No duplicate atom types in table column 0.
- Same cancel/reset/error pattern as the Descriptors run button.

---

## 6. Worker Threads

All workers live in `workers.py` and inherit from `QThread`.

### 6.1 DescriptorsWorker

```python
class DescriptorsWorker(QThread):
    progress = Signal(int)          # 0–100
    finished = Signal()
    error = Signal(str)

    def __init__(self, params: dict):
        super().__init__()
        self._params = params
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        try:
            # call qdescp_xtb_workflow with self._params
            # emit self.progress(n) at meaningful checkpoints
            # check self._cancelled between steps
            self.finished.emit()
        except Exception as e:
            self.error.emit(str(e))
```

### 6.2 NMRWorker

Identical structure, calls `qdescp_nmr_workflow`.

### 6.3 Worker lifecycle in a tab

```python
def _on_run_clicked(self):
    if self._worker and self._worker.isRunning():
        self._worker.cancel()
        self._worker.wait()
        self._reset_run_button()
        return

    params = self._collect_params()
    if not self._validate(params):
        return

    self._worker = DescriptorsWorker(params)
    self._worker.progress.connect(self._progress_bar.setValue)
    self._worker.finished.connect(self._on_finished)
    self._worker.error.connect(self._on_error)
    self._run_btn.setText("Cancel")
    self._progress_bar.setValue(0)
    self._worker.start()
```

> **Agent instruction:** Always call `worker.wait()` before discarding the reference. Dangling running threads cause segfaults on window close. Connect `QDESCPWindow.closeEvent` to cancel and wait on both workers.

---

## 7. Implementation Checklist for the Coding Agent

Work through this list **in order**. Do not skip ahead.

### Phase 1 — Window Shell

- [ ] Create `qdescp_window.py` with `QDESCPWindow(QMainWindow)`.
- [ ] Window opens via `.show()` (non-blocking). Verify the hub remains interactive while `QDESCPWindow` is open.
- [ ] Guard against duplicate windows in the hub's `open_qdescp` slot using `isVisible()`.
- [ ] `QDESCPWindow` holds a `QTabWidget` as its central widget. No other domain widgets live here.
- [ ] `QStatusBar` is created and set. `set_status(msg: str)` method is implemented.
- [ ] `closeEvent` is overridden to cancel and `.wait()` on any running workers from both tabs before accepting the close event.

### Phase 2 — Descriptors Tab Skeleton

- [ ] Create `descriptors_tab.py` with `DescriptorsTab(QWidget)`.
- [ ] Constructor accepts `status_callback: Callable[[str], None]` and stores it.
- [ ] File list (`QListWidget`) renders with Add/Remove buttons. Drag-and-drop works for `.sdf/.xyz/.pdb/.csv`.
- [ ] Mode detection (`_detect_mode`) runs on every file list change and updates the info label correctly.
- [ ] All Settings widgets render with correct defaults (verify visually).
- [ ] Advanced `QGroupBox` is checkable, starts collapsed, and shows/hides its children correctly.
- [ ] Progress bar and Run button are outside the scroll area, pinned to the bottom.

### Phase 3 — Descriptors Tab Logic

- [ ] `_collect_params()` reads all widgets and returns a `dict` matching `qdescp_xtb_workflow`'s signature.
- [ ] `_validate(params)` checks: file list non-empty; `qdescp_atoms` entries are non-empty strings if provided. Shows `QMessageBox.warning` on failure and returns `False`.
- [ ] `DescriptorsWorker` in `workers.py` is implemented. `cancel()` sets a flag; `run()` checks it between steps.
- [ ] Run button starts worker, toggles to `"Cancel"`, disables parameter widgets.
- [ ] Clicking `"Cancel"` while running calls `worker.cancel()`, waits, resets UI.
- [ ] `finished` signal resets button, sets progress to 100%, calls `status_callback`.
- [ ] `error` signal resets button, shows `QMessageBox.critical` with the error string.

### Phase 4 — NMR Tab Skeleton

- [ ] Create `nmr_tab.py` with `NMRTab(QWidget)`. Same constructor signature as `DescriptorsTab`.
- [ ] File list renders, filtered to `.sdf/.xyz` only.
- [ ] Settings widgets render with correct defaults.
- [ ] `QTableWidget` renders with three columns: Atom type, Slope, Intercept.
- [ ] Column 0 uses a `QComboBox` delegate (atom symbols). Columns 1 & 2 use `QDoubleSpinBox` delegates (4 d.p.).
- [ ] Add/Remove row buttons work. Table pre-populates with one row for the selected nucleus.
- [ ] Advanced `QGroupBox` is collapsible and starts collapsed.

### Phase 5 — NMR Tab Logic

- [ ] `_collect_params()` extracts `nmr_atoms`, `nmr_slope`, `nmr_intercept` as parallel lists from the table.
- [ ] `_validate(params)` checks: file list non-empty; table has ≥1 row; no duplicate atom types. Shows `QMessageBox.warning` on failure.
- [ ] `NMRWorker` implemented identically to `DescriptorsWorker` but calls `qdescp_nmr_workflow`.
- [ ] Run button lifecycle is identical to Descriptors (start → Cancel toggle → reset on finish/error).

### Phase 6 — Integration & Edge Cases

- [ ] Both tabs are instantiated in `QDESCPWindow` and added to the `QTabWidget` with labels `"Descriptors"` and `"NMR"`.
- [ ] Switching tabs while a worker is running on the other tab does not cancel or reset that worker.
- [ ] Closing `QDESCPWindow` while a worker is running: `closeEvent` cancels the worker, calls `.wait()`, then accepts the event. The hub is unaffected.
- [ ] Opening `QDESCPWindow` a second time from the hub while it is already open: raises and activates the existing window instead of creating a new one.
- [ ] Re-opening after closing creates a fresh instance with all widgets in their default state.
- [ ] `status_callback` correctly pushes messages to the `QStatusBar` from whichever tab calls it.

### Phase 7 — Polish

- [ ] `qdescp_atoms` `QToolButton` shows the help message box correctly.
- [ ] Solvent `QComboBox` allows freetext entry (`setEditable(True)`).
- [ ] All `QDoubleSpinBox` / `QSpinBox` widgets have correct ranges — no widget should accept logically impossible values (e.g., negative temperature, zero processes).
- [ ] Tab order (`setTabOrder`) moves focus through widgets in a logical top-to-bottom sequence within each tab.
- [ ] Minimum window size (`setMinimumSize`) prevents widgets from overlapping when the window is resized small.
- [ ] No layout warnings appear in the console on startup (check for `QLayout: Attempting to add QLayout "" to QWidget ""` style errors).

