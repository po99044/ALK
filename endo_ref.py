
def exon(a):
    
    dic ={
    "ALK_E20" : "accacccacctgcagtgtaccgccggaagcaccaggagctgcaagccatgcagatggagctgcagagccctgagtacaagctgagcaagctccgcacctcgaccatcatgaccgactacaaccccaactactgctttgctggcaagacctcctccatcagtgacctgaaggaggtgccgcggaaaaacatcaccctcattcggtgagcgccc".upper(),

    "ALK_E21" : "ttctcctttgcacaggggtctgggccatggcgcctttggggaggtgtatgaaggccaggtgtccggaatgcccaacgacccaagccccctgcaagtggctgtgaaggtaagaagtg".upper(),

    "ALK_E22" : "cccttctctgcccagacgctgcctgaagtgtgctctgaacaggacgaactggatttcctcatggaagccctgatcatcaggtaaagccac".upper(),

    "ALK_E23" : "ctctctgctctgcagcaaattcaaccaccagaacattgttcgctgcattggggtgagcctgcaatccctgccccggttcatcctgctggagctcatggcggggggagacctcaagtccttcctccgagagacccgccctcgcccggtgagtgaga".upper(),

    "ALK_E24" : "tctgtctccccacagagccagccctcctccctggccatgctggaccttctgcacgtggctcgggacattgcctgtggctgtcagtatttggaggaaaaccacttcatccaccggtgagtcaaa".upper(),

    "ALK_E25" : "ttcctttcttcccagagacattgctgccagaaactgcctcttgacctgtccaggccctggaagagtggccaagattggagacttcgggatggcccgagacatctacaggtgagtaaag".upper(),

    "ALK_E26" : "tctccttccccacagggcgagctactatagaaagggaggctgtgccatgctgccagttaagtggatgcccccagaggccttcatggaaggaatattcacttctaaaacagacacatggtaagtcagc".upper(),
    "ALK_E27" : "ctgtcccatgcccaggtcctttggagtgctgctatgggaaatcttttctcttggatatatgccataccccagcaaaagcaaccaggaagttctggagtttgtcaccagtggaggccggatggacccacccaagaactgccctgggcctgtgtatgactct".upper(),
    "ALK_E28" : "tgcttcttcttttagataccggataatgactcagtgctggcaacatcagcctgaagacaggcccaactttgccatcattttggagaggattgaatactgcacccaggtaaaacatt".upper()


    
    }

    return dic["ALK_E"+str(a)]
