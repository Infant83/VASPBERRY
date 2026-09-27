"""Exact-source metadata restoration must never alter numerical charge records."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SCRIPT=ROOT/'examples/materials/graphene-spin-chern/prepare_vasp.py'
spec=importlib.util.spec_from_file_location('graphene_prepare',SCRIPT)
prepare=importlib.util.module_from_spec(spec);spec.loader.exec_module(prepare)

class GrapheneChargeHeader(unittest.TestCase):
    def test_preserves_all_density_and_augmentation_bytes(self):
        pos=SCRIPT.parent/'inputs/POSCAR'
        lines=pos.read_text().splitlines()
        for i in (2,3,4,8,9):
            lines[i]=' '.join(f'{float(x):.6f}' for x in lines[i].split())
        payload=b' \n 36 36 300\n .12345678901E-04\naugmentation occupancies   1  12\n .33333333333E+00\n'
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'CHGCAR';out=Path(temp)/'repaired'
            source.write_bytes(('\n'.join(lines)+'\n').encode()+payload)
            before=source.read_bytes()
            record=prepare.restore_exact_charge_header(source,out,pos)
            self.assertEqual(source.read_bytes(),before)
            self.assertEqual(out.read_bytes(),pos.read_bytes()+payload)
            self.assertNotEqual(record['original_sha256'],record['output_sha256'])
    def test_rejects_different_geometry(self):
        pos=SCRIPT.parent/'inputs/POSCAR'
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'CHGCAR';out=Path(temp)/'repaired'
            source.write_text(pos.read_text().replace('0.3333333333333333','0.3334333333333333'))
            with self.assertRaisesRegex(ValueError,'geometry mismatch'):
                prepare.restore_exact_charge_header(source,out,pos)
            self.assertFalse(out.exists())
    def test_rejects_nan_geometry_and_incomplete_header(self):
        pos=SCRIPT.parent/'inputs/POSCAR'
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'CHGCAR';out=Path(temp)/'repaired'
            source.write_text(pos.read_text().replace('0.3333333333333333','nan'))
            with self.assertRaisesRegex(ValueError,'Nonfinite'):
                prepare.restore_exact_charge_header(source,out,pos)
            source.write_text('truncated\n')
            with self.assertRaisesRegex(ValueError,'complete'):
                prepare.restore_exact_charge_header(source,out,pos)
            self.assertFalse(out.exists())
    def test_rejects_close_geometry_that_is_not_six_decimal_rounding(self):
        pos=SCRIPT.parent/'inputs/POSCAR'
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)/'CHGCAR';out=Path(temp)/'repaired'
            source.write_text(pos.read_text().replace('0.3333333333333333','0.333334'))
            with self.assertRaisesRegex(ValueError,'geometry mismatch'):
                prepare.restore_exact_charge_header(source,out,pos)
            self.assertFalse(out.exists())
    def test_source_manifest_rejects_changed_inputs_and_requires_record(self):
        with tempfile.TemporaryDirectory() as temp:
            source=Path(temp)
            for name,text in [('POSCAR','structure'),('POTCAR','licensed fixture identity only'),
                              ('INCAR','ICHARG = 2\n'),('KPOINTS','mesh'),('CHGCAR','charge'),
                              ('OUTCAR','EDIFF is reached\nGeneral timing and accounting\n')]:
                (source/name).write_text(text)
            potential=prepare.sha(source/'POTCAR');structure=prepare.sha(source/'POSCAR')
            with self.assertRaisesRegex(ValueError,'requires input_manifest'):
                prepare.validate_scf_source(source,potential,structure)
            manifest={'stage':'scf','inputs':{n:prepare.sha(source/n) for n in ['POSCAR','POTCAR','INCAR','KPOINTS']}}
            (source/'input_manifest.json').write_text(json.dumps(manifest))
            record=prepare.validate_scf_source(source,potential,structure)
            self.assertEqual(record['files_sha256']['CHGCAR'],prepare.sha(source/'CHGCAR'))
            (source/'KPOINTS').write_text('different mesh')
            with self.assertRaisesRegex(ValueError,'KPOINTS'):
                prepare.validate_scf_source(source,potential,structure)

if __name__=='__main__':unittest.main()
