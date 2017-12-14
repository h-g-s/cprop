
maximize
    obj: x1 + 1.1 x2 + 1.2 x3 + x4 + x5 + x6

subject to
    r1: -x1 -x2 +x3 >= -1.0
    r2: -x3 + x4 >= 0
    r3: -x3 + x5 >= 0
    x4: -x4 + x6 >= 0
    r5: -x5 -x6 >= -1.0

binaries
    x1
    x2
    x3
    x4
    x5
    x6

End

