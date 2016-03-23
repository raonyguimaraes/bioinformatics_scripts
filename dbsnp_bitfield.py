def getBitsFromVP2(infoarr):
    for item in infoarr:
        if item.startswith("VP="):
            infoParts = item[3:]
            F0 = infoParts[0:2]
            F1_1 = infoParts[2:4]
            F1_2 = infoParts[4:6]
            F2_1 = infoParts[6:8]
            F2_2 = infoParts[8:10]
            F3 = infoParts[10:12]
            F4 = infoParts[12:14]
            F5 = infoParts[14:16]
            F6 = infoParts[16:18]
            F7 = infoParts[18:20]
            F8 = infoParts[20:22]
            F9 = infoParts[22:24]
            return "".join([bin(int(F0,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F1_1,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F1_2,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F2_1,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F2_2,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F3,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F4,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F5,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F6,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F7,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F8,16))[2:].rjust(8,"0")[::-1],
                            bin(int(F9,16))[2:].rjust(8,"0")[::-1]])
    return ""

vpBits = getBitsFromVP2(["VP=050060000a01000002110100"])

bytesArray = [vpBits[i:i+8] for i in range(0,96,8)]
print bytesArray