"""Worker for running CMIN in background thread"""
import os
import re
import shutil
import select
import sys
import tempfile
import threading
import traceback

from PySide6.QtCore import QObject, Signal
from aqme.cmin import cmin


class StreamCapture(threading.Thread):
    """Capture text from a file descriptor and forward it to Qt."""

    def __init__(self, fd, signal, stop_event):
        super().__init__(daemon=True)
        self.fd = fd
        self.signal = signal
        self.stop_event = stop_event

    def run(self):
        buffer = ""
        while not self.stop_event.is_set():
            try:
                if not select.select([self.fd], [], [], 0.1)[0]:
                    continue

                data = os.read(self.fd, 4096)
                if not data:
                    break

                text = data.decode('utf-8', errors='replace').replace('\r', '\n')
                buffer += text
                lines = buffer.split('\n')
                buffer = lines.pop()

                for line in lines:
                    line = line.rstrip()
                    if line:
                        self.signal.emit(line)
            except Exception:
                break

        tail = buffer.rstrip()
        if tail:
            self.signal.emit(tail)


class CMINWorker(QObject):
    """Worker to run CMIN optimization in background"""

    ANI_SUPPORTED_ATOMIC_NUMBERS = {
        'ani1ccx': {1, 6, 7, 8},
        'ani1x': {1, 6, 7, 8},
        'ani2x': {1, 6, 7, 8, 9, 16, 17},
    }
    ATOMIC_NUMBER_TO_SYMBOL = {
        1: 'H',
        2: 'He',
        3: 'Li',
        4: 'Be',
        5: 'B',
        6: 'C',
        7: 'N',
        8: 'O',
        9: 'F',
        10: 'Ne',
        11: 'Na',
        12: 'Mg',
        13: 'Al',
        14: 'Si',
        15: 'P',
        16: 'S',
        17: 'Cl',
        18: 'Ar',
        19: 'K',
        20: 'Ca',
        35: 'Br',
        53: 'I',
    }
    ELEMENT_SYMBOL_TO_ATOMIC_NUMBER = {
        symbol: atomic_number for atomic_number, symbol in ATOMIC_NUMBER_TO_SYMBOL.items()
    }

    ENERGY_PROPERTY_KEYS = (
        "Energy",
        "energy",
        "E_tot",
        "E",
        "SCF_Energy",
        "SCF energy",
        "Total Energy",
    )
    
    # Signals
    started = Signal()
    progress = Signal(str)  # Progress message
    result = Signal(str)  # Live output line
    finished = Signal(dict)  # Results dictionary
    error = Signal(str)  # Error message
    
    def __init__(self, files, parameters):
        super().__init__()
        self.files = files
        self.parameters = parameters
        self._is_running = False
        
    def run(self):
        """Execute CMIN with provided parameters"""
        self._is_running = True
        self.started.emit()
        staging_dir = None
        stdout_pipe_read = None
        stdout_pipe_write = None
        stderr_pipe_read = None
        stderr_pipe_write = None
        original_stdout_fd = None
        original_stderr_fd = None
        capture_stop_event = threading.Event()
        stdout_capture = None
        stderr_capture = None
        
        try:
            normalized_files = self._normalize_files(self.files)
            self.progress.emit(f"Starting CMIN optimization on {len(normalized_files)} file(s)...")

            if not normalized_files:
                raise ValueError("No valid input files were provided.")

            # AQME rewrites input paths relative to os.getcwd(); use a stable common
            # root and relative file list to prevent malformed path concatenations.
            aqme_files, run_dir, staging_dir = self._prepare_aqme_paths(normalized_files)
            output_dir = self.parameters.get('w_dir_main', '')
            output_dir = os.path.abspath(output_dir) if output_dir else ''
            
            # Build CMIN arguments
            cmin_params = dict(self.parameters)
            cmin_params.pop('w_dir_main', None)
            cmin_kwargs = {
                'files': aqme_files,
                **cmin_params,
            }

            # destination is the AQME argument for custom output location
            if output_dir:
                cmin_kwargs['destination'] = output_dir

            if str(cmin_kwargs.get('program', '')).lower() == 'ani':
                if self._ensure_torchani_compatibility():
                    self.progress.emit("Applied TorchANI compatibility patch for ANI.")

            self._validate_backend_dependencies(str(cmin_kwargs.get('program', '')).lower())
            self._validate_ani_input_elements(
                str(cmin_kwargs.get('program', '')).lower(),
                normalized_files,
                cmin_kwargs.get('ani_method'),
            )

            stdout_pipe_read, stdout_pipe_write = os.pipe()
            stderr_pipe_read, stderr_pipe_write = os.pipe()
            original_stdout_fd = os.dup(sys.stdout.fileno())
            original_stderr_fd = os.dup(sys.stderr.fileno())

            stdout_capture = StreamCapture(stdout_pipe_read, self.result, capture_stop_event)
            stderr_capture = StreamCapture(stderr_pipe_read, self.result, capture_stop_event)
            stdout_capture.start()
            stderr_capture.start()

            # Run CMIN
            previous_cwd = os.getcwd()
            try:
                os.dup2(stdout_pipe_write, sys.stdout.fileno())
                os.dup2(stderr_pipe_write, sys.stderr.fileno())
                os.chdir(run_dir)
                _ = cmin(**cmin_kwargs)
            finally:
                os.chdir(previous_cwd)
                if original_stdout_fd is not None:
                    os.dup2(original_stdout_fd, sys.stdout.fileno())
                if original_stderr_fd is not None:
                    os.dup2(original_stderr_fd, sys.stderr.fileno())
                for fd in (stdout_pipe_write, stderr_pipe_write):
                    if fd is not None:
                        try:
                            os.close(fd)
                        except OSError:
                            pass
                if stdout_capture is not None:
                    stdout_capture.join(timeout=1)
                if stderr_capture is not None:
                    stderr_capture.join(timeout=1)
                capture_stop_event.set()
                for fd in (
                    original_stdout_fd,
                    original_stderr_fd,
                    stdout_pipe_read,
                    stderr_pipe_read,
                ):
                    if fd is not None:
                        try:
                            os.close(fd)
                        except OSError:
                            pass
            
            # Collect results
            self.progress.emit("CMIN calculation finished; collecting results...")
            results = self._collect_results(output_dir, run_dir)
            
            self.progress.emit("CMIN optimization complete")
            self.finished.emit(results)
            
        except Exception as e:
            error_msg = f"CMIN failed: {str(e)}\n{traceback.format_exc()}"
            self.error.emit(error_msg)
        finally:
            if staging_dir:
                shutil.rmtree(staging_dir, ignore_errors=True)
            self._is_running = False
    
    def _collect_results(self, output_dir, run_dir):
        """Collect CMIN output statistics and energies using AQME output conventions."""
        normalized_inputs = self._normalize_files(self.files)
        base_dir = output_dir if output_dir else run_dir
        run_program = str(self.parameters.get('program', '')).strip().lower()
        discovered_outputs = self._discover_output_files(output_dir, run_dir)
        input_paths_set = {os.path.abspath(path) for path in normalized_inputs}
        discovered_outputs = [
            path for path in discovered_outputs if os.path.abspath(path) not in input_paths_set
        ]
        outputs_by_input_key, unmatched_outputs = self._classify_cmin_outputs(discovered_outputs, run_program)

        input_records = []
        for input_file in normalized_inputs:
            input_key = self._normalize_stem(os.path.splitext(os.path.basename(input_file))[0])
            output_info = outputs_by_input_key.get(input_key, {})

            filtered_file = output_info.get('filtered_file')
            filtered_count = output_info.get('filtered_count', 0)
            filtered_energies = sorted(output_info.get('filtered_energies', []))
            all_confs_count = output_info.get('all_confs_count')
            execution_status = 'success' if (filtered_file or output_info.get('all_confs_file')) else 'failed'

            input_conformers = self._count_conformers(input_file)
            if all_confs_count is None and execution_status == 'success':
                all_confs_count = input_conformers

            if all_confs_count is not None and execution_status == 'success':
                eliminated_conformers = max(0, all_confs_count - filtered_count)
            else:
                eliminated_conformers = None

            if execution_status == 'failed':
                status = 'failed'
            elif filtered_count <= 0:
                status = 'eliminated'
            elif all_confs_count is not None and filtered_count < all_confs_count:
                status = 'partial'
            else:
                status = 'success'

            input_records.append({
                'input_file': input_file,
                'input_name': os.path.basename(input_file),
                'match_key': input_key,
                'input_conformers': input_conformers,
                'all_confs_count': all_confs_count,
                'matched_output_files': [filtered_file] if filtered_file else [],
                'all_confs_file': output_info.get('all_confs_file'),
                'output_conformers': filtered_count if execution_status == 'success' else None,
                'eliminated_conformers': eliminated_conformers,
                'energies': filtered_energies,
                'min_energy': filtered_energies[0] if filtered_energies else None,
                'max_energy': filtered_energies[-1] if filtered_energies else None,
                'structures': output_info.get('filtered_structures', []),
                'execution_status': execution_status,
                'filter_outcome': self._determine_filter_outcome(filtered_count, all_confs_count, execution_status),
                'status': status,
            })

        matched_output_files = []
        for input_record in input_records:
            matched_output_files.extend(input_record['matched_output_files'])

        energy_points = []
        for input_record in input_records:
            file_label = input_record['input_name']
            for energy in input_record['energies']:
                energy_points.append({
                    'input_file': input_record['input_file'],
                    'input_name': file_label,
                    'energy': energy,
                })

        input_conformer_values = [record['all_confs_count'] for record in input_records if record['all_confs_count'] is not None]
        input_conformer_count = sum(input_conformer_values) if input_conformer_values else 0
        output_conformer_count = sum(
            record['output_conformers'] or 0 for record in input_records
        )
        conformer_level_eliminated = 0
        for record in input_records:
            if record['eliminated_conformers'] is not None:
                conformer_level_eliminated += record['eliminated_conformers']
        failed_file_count = sum(1 for record in input_records if record['execution_status'] == 'failed')
        file_level_eliminated = sum(
            1
            for record in input_records
            if record['execution_status'] == 'success' and (record['output_conformers'] or 0) <= 0
        )

        warnings = []
        if unmatched_outputs:
            warnings.append(
                f"{len(unmatched_outputs)} output file(s) did not match any selected input and were excluded from summary counts."
            )
        if any(record['all_confs_count'] is None for record in input_records):
            warnings.append("Some all-conformer files were not found. Conformer-level metrics may be incomplete.")

        results = {
            'output_files': sorted(set(matched_output_files)),
            'input_count': len(input_records),
            'output_count': len(set(matched_output_files)),
            'eliminated_count': file_level_eliminated,
            'output_dir': base_dir,
            'file_level_eliminated': file_level_eliminated,
            'failed_file_count': failed_file_count,
            'conformer_level_eliminated': conformer_level_eliminated,
            'input_conformer_count': input_conformer_count,
            'output_conformer_count': output_conformer_count,
            'per_file_data': input_records,
            'energy_values': sorted(point['energy'] for point in energy_points),
            'energy_points': energy_points,
            'unmatched_output_files': unmatched_outputs,
            'warnings': warnings,
        }

        return results

    def _classify_cmin_outputs(self, output_files, run_program):
        """Classify AQME CMIN outputs into filtered/all-conformer records per input key."""
        outputs_by_input_key = {}
        unmatched_outputs = []

        for output_file in output_files:
            if os.path.splitext(output_file)[1].lower() != '.sdf':
                unmatched_outputs.append(output_file)
                continue

            parse_info = self._parse_cmin_output_name(output_file)
            if parse_info is None:
                unmatched_outputs.append(output_file)
                continue

            if run_program and parse_info['program'] != run_program:
                continue

            input_key = parse_info['input_key']
            bucket = outputs_by_input_key.setdefault(
                input_key,
                {
                    'filtered_file': None,
                    'filtered_count': 0,
                    'filtered_energies': [],
                    'filtered_structures': [],
                    'all_confs_file': None,
                    'all_confs_count': None,
                },
            )

            include_structures = parse_info['kind'] == 'filtered'
            output_count, output_energies, output_structures = self._parse_output_file(
                output_file,
                include_structures=include_structures,
            )

            if parse_info['kind'] == 'filtered':
                bucket['filtered_file'] = output_file
                bucket['filtered_count'] = output_count
                bucket['filtered_energies'] = output_energies
                bucket['filtered_structures'] = output_structures
            elif parse_info['kind'] == 'all_confs':
                bucket['all_confs_file'] = output_file
                bucket['all_confs_count'] = output_count

        return outputs_by_input_key, unmatched_outputs

    def _parse_cmin_output_name(self, file_path):
        """Parse AQME CMIN output name and return kind/program/input mapping info."""
        stem = os.path.splitext(os.path.basename(file_path))[0]

        all_match = re.match(r'^(?P<input>.+)_(?P<program>xtb|ani)_all_confs$', stem, flags=re.IGNORECASE)
        if all_match:
            return {
                'kind': 'all_confs',
                'program': all_match.group('program').lower(),
                'input_key': self._normalize_stem(all_match.group('input')),
            }

        filtered_match = re.match(r'^(?P<input>.+)_(?P<program>xtb|ani)$', stem, flags=re.IGNORECASE)
        if filtered_match:
            return {
                'kind': 'filtered',
                'program': filtered_match.group('program').lower(),
                'input_key': self._normalize_stem(filtered_match.group('input')),
            }

        return None

    def _normalize_stem(self, text):
        """Normalize a filename stem to a stable comparison key."""
        return re.sub(r'[^a-z0-9]+', '', str(text).lower())

    def _discover_output_files(self, output_dir, run_dir):
        """Find output files in common AQME layouts without duplicates."""
        discovered = []
        search_dirs = []
        if output_dir:
            search_dirs.extend([output_dir, os.path.join(output_dir, 'CMIN')])
        search_dirs.extend([run_dir, os.path.join(run_dir, 'CMIN')])

        visited_dirs = set()
        for cmin_dir in search_dirs:
            if not cmin_dir:
                continue
            abs_dir = os.path.abspath(cmin_dir)
            if abs_dir in visited_dirs or not os.path.exists(abs_dir):
                continue
            visited_dirs.add(abs_dir)

            for root, _, files in os.walk(abs_dir):
                for file in files:
                    if file.endswith('.sdf') or file.endswith('.xyz'):
                        discovered.append(os.path.abspath(os.path.join(root, file)))

        return sorted(set(discovered))

    def _make_match_key(self, file_path):
        """Create a normalized key used to map input and output files."""
        stem = os.path.splitext(os.path.basename(file_path))[0]
        return self._normalize_stem(stem)

    def _match_output_to_input(self, output_key, input_records):
        """Best-effort mapping of an output file key to one input record index."""
        if not output_key:
            return None

        best_index = None
        best_score = -1
        for index, record in enumerate(input_records):
            input_key = record.get('match_key', '')
            if not input_key:
                continue
            if output_key == input_key:
                score = len(input_key) + 1000
            elif output_key.startswith(input_key) or input_key in output_key:
                score = len(input_key)
            elif input_key.startswith(output_key):
                score = len(output_key)
            else:
                continue

            if score > best_score:
                best_score = score
                best_index = index

        return best_index

    def _parse_output_file(self, file_path, include_structures=False):
        """Return output conformer count and parsed energies for one output file."""
        ext = os.path.splitext(file_path)[1].lower()
        if ext == '.sdf':
            count, energies, structures = self._parse_sdf_summary(file_path)
            if include_structures and count > 0:
                enriched = self._enrich_structures_with_rdkit(file_path)
                if enriched:
                    if not energies:
                        energies = [
                            float(item['energy'])
                            for item in enriched
                            if item.get('energy') is not None
                        ]
                    return count, energies, enriched
            return count, energies, structures
        if ext == '.xyz':
            count = self._count_xyz_conformers(file_path)
            return count, [], []
        return 0, [], []

    def _parse_sdf_file(self, file_path):
        """Read conformer count and energies from an SDF file."""
        count, energies, _ = self._parse_sdf_summary(file_path)
        return count, energies

    def _parse_sdf_summary(self, file_path):
        """Read conformer count, energies, and lightweight structure rows from one SDF file.

        This parser intentionally avoids RDKit conformer iteration and RMSD calculations,
        which keeps result collection fast for large SDF outputs.
        """
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as handle:
                content = handle.read()
        except Exception:
            return 0, [], []

        blocks = [block.strip() for block in content.split('$$$$')]
        entries = []
        for index, block in enumerate(blocks):
            if not block:
                continue

            lines = [line.rstrip('\n') for line in block.splitlines()]
            name = lines[0].strip() if lines and lines[0].strip() else f"Conformer {index + 1}"
            energy = self._extract_energy_from_sdf_block(lines)
            entries.append(
                {
                    'name': name,
                    'energy': energy,
                }
            )

        energies = [float(entry['energy']) for entry in entries if entry['energy'] is not None]
        sortable_entries = [entry for entry in entries if entry['energy'] is not None]
        unsorted_entries = [entry for entry in entries if entry['energy'] is None]
        sortable_entries.sort(key=lambda item: item['energy'])
        ordered_entries = sortable_entries + unsorted_entries

        structures = []
        for rank, entry in enumerate(ordered_entries, start=1):
            structures.append(
                {
                    'rank': rank,
                    'name': entry['name'],
                    'energy': entry['energy'],
                    'rmsd_to_best': 0.0 if rank == 1 else None,
                }
            )

        return len(entries), energies, structures

    def _extract_energy_from_sdf_block(self, lines):
        """Extract an energy value from one SDF text block."""
        energy_keys = {str(key).strip().lower() for key in self.ENERGY_PROPERTY_KEYS}
        index = 0
        total_lines = len(lines)
        while index < total_lines:
            line = lines[index].strip()
            match = re.match(r'^>\s*<([^>]+)>\s*$', line)
            if match:
                property_name = match.group(1).strip().lower()
                value_index = index + 1
                while value_index < total_lines and not lines[value_index].strip():
                    value_index += 1
                if (property_name in energy_keys or 'energy' in property_name) and value_index < total_lines:
                    try:
                        return float(str(lines[value_index]).strip())
                    except (TypeError, ValueError):
                        pass
            index += 1

        return None

    def _enrich_structures_with_rdkit(self, file_path):
        """Build detailed structure rows (including RMSD) for filtered outputs only."""
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolAlign
        except Exception:
            return []

        mol_entries = []
        try:
            supplier = Chem.SDMolSupplier(file_path, removeHs=False, sanitize=False)
            for index, mol in enumerate(supplier):
                if mol is None:
                    continue
                energy = self._extract_energy_from_mol(mol)
                name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"Conformer {index + 1}"
                mol_entries.append(
                    {
                        'name': str(name).strip() or f"Conformer {index + 1}",
                        'energy': energy,
                        'mol': mol,
                    }
                )
        except Exception:
            return []

        if not mol_entries:
            return []

        sortable = [entry for entry in mol_entries if entry['energy'] is not None]
        unsorted = [entry for entry in mol_entries if entry['energy'] is None]
        sortable.sort(key=lambda item: item['energy'])
        ordered_entries = sortable + unsorted

        reference_mol = sortable[0]['mol'] if sortable else ordered_entries[0]['mol']
        structures = []
        for rank, entry in enumerate(ordered_entries, start=1):
            rmsd_value = None
            try:
                rmsd_value = float(rdMolAlign.GetBestRMS(reference_mol, entry['mol']))
            except Exception:
                rmsd_value = 0.0 if entry['mol'] is reference_mol else None

            structures.append(
                {
                    'rank': rank,
                    'name': entry['name'],
                    'energy': entry['energy'],
                    'rmsd_to_best': rmsd_value,
                }
            )

        return structures

    def _determine_filter_outcome(self, output_conformers, all_confs_count, execution_status):
        """Return filtering outcome label independent from execution status."""
        if execution_status != 'success':
            return 'execution_failed'
        if output_conformers <= 0:
            return 'eliminated'
        if all_confs_count is None:
            return 'retained_unknown'
        if output_conformers < all_confs_count:
            return 'filtered'
        return 'retained_all'

    def _extract_energy_from_mol(self, mol):
        """Extract one numeric energy value from known SDF properties."""
        for key in self.ENERGY_PROPERTY_KEYS:
            if not mol.HasProp(key):
                continue
            value = mol.GetProp(key)
            try:
                return float(str(value).strip())
            except (TypeError, ValueError):
                continue
        return None

    def _count_conformers(self, file_path):
        """Count conformers in input/output files using extension-aware parsing."""
        ext = os.path.splitext(file_path)[1].lower()
        if ext == '.sdf':
            count, _ = self._parse_sdf_file(file_path)
            return count
        if ext == '.xyz':
            return self._count_xyz_conformers(file_path)
        return 1

    def _count_xyz_conformers(self, file_path):
        """Count XYZ blocks in one file by iterating atom-count headers."""
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as handle:
                lines = handle.readlines()
        except Exception:
            return None

        conformers = 0
        index = 0
        total_lines = len(lines)
        while index < total_lines:
            stripped = lines[index].strip()
            if not stripped:
                index += 1
                continue
            try:
                atom_count = int(stripped)
            except ValueError:
                index += 1
                continue

            block_end = index + atom_count + 2
            if atom_count > 0 and block_end <= total_lines:
                conformers += 1
                index = block_end
            else:
                index += 1

        return conformers

    def _normalize_files(self, files):
        """Return unique, existing absolute file paths preserving order."""
        normalized = []
        seen = set()
        for file in files or []:
            if not file:
                continue
            abs_path = os.path.abspath(str(file))
            if abs_path in seen:
                continue
            if os.path.isfile(abs_path):
                seen.add(abs_path)
                normalized.append(abs_path)
        return normalized

    def _prepare_aqme_paths(self, normalized_files):
        """Build AQME-compatible file arguments and execution directory.

        AQME path and conversion helpers are more stable when files are flat
        basenames in the current working directory.
        """
        if not normalized_files:
            return [], os.getcwd(), None

        basenames = [os.path.basename(path) for path in normalized_files]
        duplicate_names = sorted({name for name in basenames if basenames.count(name) > 1})
        if duplicate_names:
            duplicate_list = ", ".join(duplicate_names)
            raise ValueError(
                "Selected files contain duplicate names from different folders. "
                f"Rename one of these files and try again: {duplicate_list}"
            )

        staging_dir = tempfile.mkdtemp(prefix="aqmeasy_cmin_")
        aqme_files = []
        output_names = set()
        generated_sdf_names = set()
        for file_path in normalized_files:
            staged_name = os.path.basename(file_path)
            staged_path = os.path.join(staging_dir, staged_name)
            shutil.copy2(file_path, staged_path)

            ext = os.path.splitext(staged_name)[1].lower()
            if ext == '.xyz':
                sdf_name = f"{os.path.splitext(staged_name)[0]}.sdf"
                if sdf_name in output_names:
                    raise ValueError(
                        "Selected files generate duplicate staged SDF names. "
                        f"Rename one input and try again: {sdf_name}"
                    )
                self._convert_xyz_to_sdf(staged_path)
                aqme_files.append(sdf_name)
                output_names.add(sdf_name)
                generated_sdf_names.add(sdf_name)
            else:
                if staged_name in output_names:
                    if staged_name in generated_sdf_names:
                        raise ValueError(
                            "Selected files generate duplicate staged SDF names. "
                            f"Rename one input and try again: {staged_name}"
                        )
                    raise ValueError(
                        "Selected files contain duplicate names from different folders. "
                        f"Rename one input and try again: {staged_name}"
                    )
                aqme_files.append(staged_name)
                output_names.add(staged_name)

        return aqme_files, staging_dir, staging_dir

    def _convert_xyz_to_sdf(self, xyz_path):
        """Convert one staged XYZ file to SDF using AQME/OpenBabel helper."""
        try:
            from aqme.csearch.utils import xyz_2_sdf
            xyz_2_sdf(xyz_path)
        except Exception as exc:
            xyz_name = os.path.basename(xyz_path)
            raise RuntimeError(
                "Failed to convert XYZ input to SDF before CMIN run. "
                f"File: {xyz_name}. Ensure OpenBabel is installed and the XYZ is valid."
            ) from exc

    def _ensure_torchani_compatibility(self):
        """Patch TorchANI ANI models with species_to_tensor when missing.

        AQME currently calls model.species_to_tensor(elements), but newer TorchANI
        versions expose species conversion through different APIs.
        """
        try:
            import torch
            from torchani.arch import ANI as TorchANIModel
            from torchani.constants import ATOMIC_NUMBER
        except Exception:
            return False

        if hasattr(TorchANIModel, 'species_to_tensor'):
            return False

        def _species_to_tensor(self, elements):
            if isinstance(elements, str):
                # Split concatenated symbols (e.g. "CClH") into ["C", "Cl", "H"]
                symbols = re.findall(r'[A-Z][a-z]?', elements)
            elif isinstance(elements, (list, tuple)):
                symbols = [str(ele) for ele in elements]
            else:
                raise TypeError(f"Unsupported elements format: {type(elements).__name__}")

            try:
                numbers = [ATOMIC_NUMBER[symbol] for symbol in symbols]
            except KeyError as exc:
                raise KeyError(f"Unsupported element symbol for ANI: {exc.args[0]}") from exc

            return torch.tensor(numbers, dtype=torch.long)

        setattr(TorchANIModel, 'species_to_tensor', _species_to_tensor)
        return True

    def _validate_backend_dependencies(self, program):
        """Fail early with a clear message when external backend binaries are missing."""
        if program == 'xtb':
            required_bins = ['crest', 'xtb']
            missing = [binary for binary in required_bins if shutil.which(binary) is None]
            if missing:
                missing_list = ', '.join(missing)
                raise FileNotFoundError(
                    "xTB backend requires external binaries not found in PATH: "
                    f"{missing_list}. Install them and ensure PATH includes their location."
                )

    def _validate_ani_input_elements(self, program, input_files, ani_method):
        """Validate ANI model element coverage before launching AQME/TorchANI."""
        if program != 'ani':
            return

        model_key = self._normalize_ani_method(ani_method)
        supported_atomic_numbers = self.ANI_SUPPORTED_ATOMIC_NUMBERS[model_key]
        input_atomic_numbers = self._collect_input_atomic_numbers(input_files)
        unsupported_atomic_numbers = sorted(input_atomic_numbers - supported_atomic_numbers)
        if not unsupported_atomic_numbers:
            return

        unsupported_elements = self._format_atomic_numbers(unsupported_atomic_numbers)
        supported_elements = self._format_atomic_numbers(sorted(supported_atomic_numbers))
        model_label = self._format_ani_method_label(model_key)
        compatible_methods = self._find_compatible_ani_methods(input_atomic_numbers)
        if compatible_methods:
            compatibility_hint = (
                f"Compatible ANI methods for this input: {', '.join(compatible_methods)}. "
                "Choose one of them or remove unsupported elements from the input."
            )
        else:
            compatibility_hint = (
                "No ANI method in AQME supports this element set. Use xTB or remove "
                "unsupported elements from the input."
            )
        raise ValueError(
            f"Selected ANI method '{model_label}' does not support element(s): {unsupported_elements}. "
            f"Supported elements for {model_label}: {supported_elements}. "
            f"{compatibility_hint}"
        )

    def _normalize_ani_method(self, ani_method):
        """Return canonical ANI method key from UI/AQME-style values."""
        normalized = re.sub(r'[^a-z0-9]+', '', str(ani_method or 'ANI2x').lower())
        aliases = {
            'ani1ccx': 'ani1ccx',
            'ani1x': 'ani1x',
            'ani2x': 'ani2x',
        }
        canonical = aliases.get(normalized)
        if canonical is None:
            supported = ', '.join(self._format_ani_method_label(key) for key in sorted(aliases.values()))
            raise ValueError(
                f"Unknown ANI method '{ani_method}'. Supported methods are: {supported}."
            )
        return canonical

    def _format_ani_method_label(self, method_key):
        """Render a canonical ANI method key as a user-facing label."""
        if method_key == 'ani1ccx':
            return 'ANI1ccx'
        if method_key == 'ani1x':
            return 'ANI1x'
        if method_key == 'ani2x':
            return 'ANI2x'
        return str(method_key)

    def _find_compatible_ani_methods(self, input_atomic_numbers):
        """Return ANI methods that support every element in the input set."""
        compatible_methods = []
        for method_key, supported_atomic_numbers in self.ANI_SUPPORTED_ATOMIC_NUMBERS.items():
            if set(input_atomic_numbers).issubset(supported_atomic_numbers):
                compatible_methods.append(self._format_ani_method_label(method_key))
        return compatible_methods

    def _collect_input_atomic_numbers(self, input_files):
        """Collect all unique atomic numbers from selected SDF/XYZ input files."""
        atomic_numbers = set()
        for file_path in input_files or []:
            ext = os.path.splitext(str(file_path))[1].lower()
            if ext == '.sdf':
                atomic_numbers.update(self._collect_atomic_numbers_from_sdf(file_path))
            elif ext == '.xyz':
                atomic_numbers.update(self._collect_atomic_numbers_from_xyz(file_path))
        return atomic_numbers

    def _collect_atomic_numbers_from_sdf(self, file_path):
        """Collect atomic numbers from an SDF file using RDKit or text fallback."""
        try:
            from rdkit import Chem
        except Exception:
            return self._collect_atomic_numbers_from_sdf_text(file_path)

        atomic_numbers = set()
        try:
            supplier = Chem.SDMolSupplier(file_path, removeHs=False, sanitize=False)
            for mol in supplier:
                if mol is None:
                    continue
                for atom in mol.GetAtoms():
                    atomic_numbers.add(int(atom.GetAtomicNum()))
        except Exception:
            return self._collect_atomic_numbers_from_sdf_text(file_path)
        return atomic_numbers

    def _collect_atomic_numbers_from_sdf_text(self, file_path):
        """Collect atomic numbers from SDF by parsing V2000/V3000 atom lines."""
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as handle:
                lines = handle.readlines()
        except Exception:
            return set()

        atomic_numbers = set()
        index = 0
        total_lines = len(lines)
        while index < total_lines:
            line = lines[index]
            if 'V3000' in line:
                index = self._parse_v3000_atom_block(lines, index, atomic_numbers)
                continue

            if index + 3 >= total_lines:
                break

            counts_line = lines[index + 3]
            try:
                atom_count = int(counts_line[:3])
            except Exception:
                index += 1
                continue

            atom_start = index + 4
            atom_end = min(atom_start + atom_count, total_lines)
            for atom_line in lines[atom_start:atom_end]:
                self._add_atomic_number_from_symbol(self._extract_v2000_symbol(atom_line), atomic_numbers)

            index = atom_end

        return atomic_numbers

    def _parse_v3000_atom_block(self, lines, start_index, atomic_numbers):
        """Parse one V3000 atom block and add atomic numbers into the accumulator."""
        index = start_index
        total_lines = len(lines)
        while index < total_lines:
            line = lines[index].strip()
            if line.startswith('M  V30 BEGIN ATOM'):
                index += 1
                while index < total_lines:
                    atom_line = lines[index].strip()
                    if atom_line.startswith('M  V30 END ATOM'):
                        return index + 1
                    if atom_line.startswith('M  V30'):
                        parts = atom_line.split()
                        if len(parts) >= 4:
                            self._add_atomic_number_from_symbol(parts[3], atomic_numbers)
                    index += 1
                return total_lines
            if line.startswith('$$$$'):
                return index + 1
            index += 1
        return total_lines

    def _extract_v2000_symbol(self, atom_line):
        """Extract element symbol from one V2000 atom line."""
        symbol = atom_line[31:34].strip() if len(atom_line) >= 34 else ''
        if symbol:
            return symbol
        parts = atom_line.split()
        if len(parts) >= 4:
            return parts[3]
        return ''

    def _collect_atomic_numbers_from_xyz(self, file_path):
        """Collect atomic numbers from XYZ symbols across all conformer blocks."""
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as handle:
                lines = handle.readlines()
        except Exception:
            return set()

        atomic_numbers = set()
        index = 0
        total_lines = len(lines)
        while index < total_lines:
            stripped = lines[index].strip()
            if not stripped:
                index += 1
                continue
            try:
                atom_count = int(stripped)
            except ValueError:
                index += 1
                continue

            atom_start = index + 2
            atom_end = min(atom_start + atom_count, total_lines)
            for atom_line in lines[atom_start:atom_end]:
                parts = atom_line.split()
                if parts:
                    self._add_atomic_number_from_symbol(parts[0], atomic_numbers)
            index = atom_end

        return atomic_numbers

    def _add_atomic_number_from_symbol(self, symbol, atomic_numbers):
        """Translate one element symbol and add it to the atomic-number set."""
        normalized_symbol = self._normalize_symbol(symbol)
        if not normalized_symbol:
            return
        atomic_number = self.ELEMENT_SYMBOL_TO_ATOMIC_NUMBER.get(normalized_symbol)
        if atomic_number is not None:
            atomic_numbers.add(atomic_number)

    def _normalize_symbol(self, symbol):
        """Normalize free-form atom symbol text to proper element symbol casing."""
        text = str(symbol or '').strip()
        if not text:
            return ''
        letters = ''.join(ch for ch in text if ch.isalpha())
        if not letters:
            return ''
        if len(letters) == 1:
            return letters.upper()
        return letters[0].upper() + letters[1:].lower()

    def _format_atomic_numbers(self, atomic_numbers):
        """Format atomic numbers as sorted element labels with atomic numbers."""
        labels = []
        for atomic_number in sorted(atomic_numbers):
            symbol = self.ATOMIC_NUMBER_TO_SYMBOL.get(atomic_number)
            if symbol is None:
                labels.append(f"Z={atomic_number}")
            else:
                labels.append(f"{symbol}(Z={atomic_number})")
        return ', '.join(labels)
    
    def stop(self):
        """Request worker to stop (CMIN doesn't support interruption, but we can flag it)"""
        self._is_running = False