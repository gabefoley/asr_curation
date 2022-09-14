import scripts.seqcurate

def test_randstring():

    default_rand = scripts.seqcurate.randstring()
    five_rand = scripts.seqcurate.randstring(5)

    assert(len(default_rand) == 10)
    assert (len(five_rand) == 5)