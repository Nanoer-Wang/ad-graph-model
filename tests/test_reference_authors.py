"""Guard against empty Crossref consortium records in author positions."""
import unittest,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class ReferenceAuthors(unittest.TestCase):
    def test_all_first_authors_are_named(self):
        audit=json.loads((ROOT/'manuscript/reference_audit.json').read_text(encoding='utf-8'))
        for ref in audit:
            self.assertTrue(ref['first_author'].strip(),ref['number'])
        vogel=next(r for r in audit if r['doi']=='10.1038/s41591-021-01309-6')
        self.assertEqual(vogel['first_author'],'Vogel JW')
if __name__=='__main__':unittest.main()
