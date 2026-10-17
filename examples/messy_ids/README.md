# messy_ids: protein and DNA IDs that don't quite match

Ten synthetic Arabidopsis-style genes. The protein alignment uses transcript IDs such as
`AT1G01010.1`, while the DNA file uses gene IDs such as `AT1G01010`, and the DNA records are
in a different order. This is very common when proteins and CDS come from different
downloads. A tool that pairs sequences by exact ID (or by position) will pair the wrong
protein with the wrong DNA. The fix is to match IDs after removing the `.1` version
suffix; once paired correctly, every gene fits the standard genetic code (1).
