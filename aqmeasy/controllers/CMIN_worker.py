"""Worker for running CMIN in background thread"""
from PySide6.QtCore import QObject, Signal
from aqme.cmin import cmin
import traceback
import os
import re
import shutil


class CMINWorker(QObject):
    """Worker to run CMIN optimization in background"""

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
        
        try:
            normalized_files = self._normalize_files(self.files)
            self.progress.emit(f"Starting CMIN optimization on {len(normalized_files)} file(s)...")

            if not normalized_files:
                raise ValueError("No valid input files were provided.")

            # AQME currently infers/rewrites working paths from os.getcwd() and files.
            # Running from the first input file directory avoids invalid path rewriting.
            run_dir = os.path.dirname(normalized_files[0])
            output_dir = self.parameters.get('w_dir_main', '')
            output_dir = os.path.abspath(output_dir) if output_dir else ''
            
            # Build CMIN arguments
            cmin_params = dict(self.parameters)
            cmin_params.pop('w_dir_main', None)
            cmin_kwargs = {
                'files': normalized_files,
                **cmin_params,
            }

            # destination is the AQME argument for custom output location
            if output_dir:
                cmin_kwargs['destination'] = output_dir

            if str(cmin_kwargs.get('program', '')).lower() == 'ani':
                if self._ensure_torchani_compatibility():
                    self.progress.emit("Applied TorchANI compatibility patch for ANI.")

            self._validate_backend_dependencies(str(cmin_kwargs.get('program', '')).lower())
            
            # Run CMIN
            previous_cwd = os.getcwd()
            try:
                os.chdir(run_dir)
                _ = cmin(**cmin_kwargs)
            finally:
                os.chdir(previous_cwd)
            
            # Collect results
            results = self._collect_results(output_dir, run_dir)
            
            self.progress.emit("CMIN optimization complete")
            self.finished.emit(results)
            
        except Exception as e:
            error_msg = f"CMIN failed: {str(e)}\n{traceback.format_exc()}"
            self.error.emit(error_msg)
        finally:
            self._is_running = False
    
    def _collect_results(self, output_dir, run_dir):
        """Collect CMIN output statistics and energies using AQME output conventions."""
        normalized_inputs = self._normalize_files(self.files)
        base_dir = output_dir if output_dir else run_dir
        run_program = str(self.parameters.get('program', '')).strip().lower()
        discovered_outputs = self._discover_output_files(output_dir, run_dir)
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

            output_count, output_energies, output_structures = self._parse_output_file(output_file)

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

    def _parse_output_file(self, file_path):
        """Return output conformer count and parsed energies for one output file."""
        ext = os.path.splitext(file_path)[1].lower()
        if ext == '.sdf':
            count, energies = self._parse_sdf_file(file_path)
            structures = self._parse_sdf_structures(file_path)
            if not structures and energies:
                # Fallback synthetic structure rows when detailed parsing is unavailable.
                structures = [
                    {
                        'rank': index + 1,
                        'name': f'Conformer {index + 1}',
                        'energy': energy,
                        'rmsd_to_best': 0.0 if index == 0 else None,
                    }
                    for index, energy in enumerate(sorted(energies))
                ]
            return count, energies, structures
        if ext == '.xyz':
            count = self._count_xyz_conformers(file_path)
            return count, [], []
        return 0, [], []

    def _parse_sdf_file(self, file_path):
        """Read conformer count and energies from an SDF file."""
        structures = self._parse_sdf_structures(file_path)
        energies = []
        for structure in structures:
            energy = structure.get('energy')
            if energy is not None:
                energies.append(float(energy))

        if structures:
            return len(structures), energies

        try:
            from rdkit import Chem
        except Exception:
            return 0, []

        count = 0
        energies = []
        try:
            supplier = Chem.SDMolSupplier(file_path, removeHs=False, sanitize=False)
            for mol in supplier:
                if mol is None:
                    continue
                count += 1
                energy = self._extract_energy_from_mol(mol)
                if energy is not None:
                    energies.append(energy)
        except Exception:
            return 0, []

        return count, energies

    def _parse_sdf_structures(self, file_path):
        """Return per-structure records including energy and RMSD to lowest-energy conformer."""
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
                        'index': index,
                        'name': str(name).strip() or f"Conformer {index + 1}",
                        'energy': energy,
                        'mol': mol,
                    }
                )
        except Exception:
            return []

        if not mol_entries:
            return []

        sortable_entries = []
        unsorted_entries = []
        for entry in mol_entries:
            if entry['energy'] is None:
                unsorted_entries.append(entry)
            else:
                sortable_entries.append(entry)

        sortable_entries.sort(key=lambda item: item['energy'])
        ordered_entries = sortable_entries + unsorted_entries

        reference_mol = sortable_entries[0]['mol'] if sortable_entries else ordered_entries[0]['mol']
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
    
    def stop(self):
        """Request worker to stop (CMIN doesn't support interruption, but we can flag it)"""
        self._is_running = False