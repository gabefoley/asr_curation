from brendapy import BrendaParser



BRENDA_PARSER = BrendaParser()

def test_parse_proteins_for_ec():
    proteins = BRENDA_PARSER.get_proteins("1.1.1.1")




