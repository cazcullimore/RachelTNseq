import unittest

from objects import Gene, HopSite


class NormalizationTests(unittest.TestCase):
    def test_fractional_sites_survive_in_gene_and_hop_tables(self):
        gene = Gene()
        gene.create_intergenic_region(1, 100, "region", 2)
        for position, counts in [(10, [1, 1]), (20, [0, 999])]:
            hop = HopSite(position, False, "+", 2)
            hop.hops = counts
            gene.hop_list.append(hop)

        # Library totals are 1 and 1000; normalize both to the smaller total.
        rows = gene.write_hops(1, [1.0, 0.001]).splitlines()
        first_site = rows[1].split("\t")
        self.assertEqual(first_site[2:4], ["1", "1"])
        self.assertEqual([float(x) for x in first_site[4:6]], [1.0, 0.001])
        self.assertEqual(gene.hop_totals(), [1, 1000])
        gene_row = gene.write_gene().split("\t")
        self.assertEqual(gene_row[2:4], ["1", "1000"])
        for total in gene_row[4:6]:
            self.assertAlmostEqual(float(total), 1.0)
        self.assertEqual(rows[0], gene.write_gene())

    def test_igv_preserves_fractional_counts_and_strand(self):
        hop = HopSite(10, True, "-", 2)
        hop.hops = [0, 1]
        for strand, sign in [("-", -1), ("+", 1)]:
            hop.strand = strand
            row = hop.get_hops_for_igv("gene", "ref", True, [1.0, 0.001])
            self.assertEqual([float(x) for x in row.split("\t")[4:]],
                             [0.0, sign * 0.001])

    def test_raw_igv_remains_integer_counts(self):
        hop = HopSite(10, True, "-", 2)
        hop.hops = [0, 1]
        row = hop.get_hops_for_igv("gene", "ref", False, [1.0, 0.001])
        self.assertEqual(row.split("\t")[4:], ["0", "-1"])


if __name__ == "__main__":
    unittest.main()
