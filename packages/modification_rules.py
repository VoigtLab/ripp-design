

aas = 'FILWVMYCPAGTSQNEDHKR'
core_rules_r = dict([
    ['tgnb' , r"[GFL][PHL][DS][TSVH][T][IMEH][VSR][T][EKRSAVLG][T][ITS][E][NSQM][AVRFW][D][PI][DS][EAMW][YLST][FYEP][LAFQG]"],
    ['plpxy', r"[AT][VTCQSNKRP][AS][A][MLN][Y][G][VATS][V][FITEKP][P]"],
    ['paap' , r"[E][E][NSTDKRYV][AMTGP][MVSTNHD][YF][TSAVLWYG][KRSAY][GSTEHMVL][QNTHALVPG][VLATKRW][IKQ][VING][LVAKQRS][SA]"],
    ['lynd' , r"[LIMVAYWFTQSNKRHEDGP][C][SNTQYWFDEKRHLIMVAGP]"],
    ['lasf' , r"[Q][LVYS][V][GAW][RV][RV][NEL][I]$"],
    ['pals' , r"[GPNSQHRVI][C][GS][GSQCH]"],
    ['epid' , r"[GL][SE][FVKG][NQTRY][SQRG][YLN][CV][C]$"],
    ['thcok', r"[YWA][S]$"],
    ['padek', r"[HYFWALGP][YLC][D][S]$"],
    ['tevp' , r"[FWMYCAGTSQNDHKR]"]]) 

core_rules = dict([(k,[set(cr) for cr in v.strip('$').strip(']').strip('[').split('][')]) for k, v in core_rules_r.items()])

#spacing rule tuple is arranged as:
# (optimal RS-to-mod distance,
#  position of the mod in the core motif above (0 is at the beginning of the motif),
#  insertion spring constant (compression, kc),
#  deletion spring constant (stretching, ks))

spacing_rules = dict([
    ['tgnb' , (37  , 4, 140  , 40  )],
    ['plpxy', (6   , 5, 16   , 1700)],
    ['paap' , (0   , 0, 5500 , 5800)],
    ['lynd' , (8   , 1, 10   , 300 )],
    ['lasf' , (None, 7, None , None)],
    ['pals' , (None, 1, None , None)],
    ['epid' , (None, 7, None , None)],
    ['thcok', (None, 1, None , None)],
    ['padek', (None, 3, None , None)],
    ['tevp' , (0   , 0, 1000000000, 1000000000)]])

recognition_sites = dict([
    ['tgnb' , "PYIAKYVEE"],
    ['plpxy', "ELNEEELEAIAG"],
    ['paap' , "FSTLSQRISAIT"],
    ['lynd' , "LAELSEEAL"],
    ['tevp' , "ENLYFQ"]])

leaders = dict([
    ['tgnb' , ("YR","QTLQNSTNLVYDDITQISFINKEKNVKKINL")],
    ['plpxy', ("SIESAKAFYQRMTDDASFRTPFEAELSKEERQQLIKDSGYDFTAEEWQQAMTEIQAARSNE","G")],
    ['paap' , ("IK","")],
    ['lynd' , ("NKKNILPQLGQPVIRLTAGQLSSQ","GGVDAS")],
    ['tevp' , ("","")]])
