def my_hook(t):
    last_b = [0]

    def inner(b=1, bsize=1, tsize=None):
        """
            b  : int, optional
                Number of blocks just transferred [default: 1].
            bsize  : int, optional
                Size of each block (in tqdm units) [default: 1].
            tsize  : int, optional
                Total size (in tqdm units). If [default: None] remains unchanged.
            """
        if tsize is not None:
            t.total = tsize
        t.update((b - last_b[0]) * bsize)
        last_b[0] = b

    return inner


BASES = {"A": "T", "T": "A", "G": "C", "C": "G", "-": "-"}


def get_complement(seq):
    return list(sorted(BASES[x] for x in seq))