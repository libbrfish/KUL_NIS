from histmatch import scale_interval

def test_scale_interval():
    v = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    assert scale_interval(v, 2, 4, 0, 2) == (1, 0, 1, 4, 5, 6, 7, 8, 9, 10)


