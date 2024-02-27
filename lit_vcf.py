#!/usr/bin/python3
# -*- coding: utf-8 -*-

#extrait les données intéressante de tous les fichiers vcf et crée deux fichiers tsv :
#   un fichier Ptout avec tous les variants
#   un fichier Pmerged avec les variants rassemblés par position, type, échantillon et passage
#nécessite d'être dans le même répertoire que tous les fichiers VCF avec les noms originaux


#listes des passages et échantillons à parcourir :
liste_passages = ("P15", "P30", "P50", "P65")
liste_echantillons = list(range(1, 11))

#filtres à appliquer, peuvent être modifiés
seuil_freq = 0.1
seuil_profondeur = 100
#nombre de nucléotides d'écart maximum entre les variants à merger:
diffmerge = 10

#dictionnaire toutes les insertions / délétions :
DtoutINS = {}
DtoutDEL = {}
#dictionnaire contenant les variants mergés:
Dmerged = {}

#ouverture des fichier de sortie :
fichierBrut = open("Ptout.tsv", 'w')
fichierMerge = open("Pmerged.tsv", 'w')

#header des fichier de sortie:
fichierBrut.write("passage\techantillon\ttype_stress\tposition\tseq_reference\tseq_alternative\ttype_sv\tfrequence\tprofondeur\n")
fichierMerge.write("passage\techantillon\ttype_stress\tposition\ttype_sv\tfrequence\tprofondeur\n")




#écrit un variant dans le fichier donné, ne retourne rien
def ecritVariant (v, fichier):  #v dictionnaire d'un variant, fichier un fichier ouvert
    for _, item in v.items():
        fichier.write(f'{item}\t')
    fichier.write('\n')

#fonction de remplissage de dictionnaire, ne retourne rien (modifie d)
def remplitDict (d, v): #d dictionnaire de tous les passages/échantillons, v dictionnaire d'un variant
    if v["position"] not in list(d[nom_passage][num_echantillon]):
        d[nom_passage][num_echantillon][v["position"]] = []
    d[nom_passage][num_echantillon][v["position"]].append(v)

#crée un variant "moyen" à partir d'une liste de variants, retourne le dictionnaire du variant 
def mergeVar (listVar, posm, indel) :   #listVar liste de dictionnaires de variants, posm entier (position moyenne des variants), indel le type de variant
    freqtot = 0
    proftot = 0
    if num_echantillon < 6 :
        type_stress = 'froid'
    else :
        type_stress = 'chaud'
    varmerged = {'passage': nom_passage, 'echantillon': num_echantillon, 'type_stress': type_stress}
    varmerged['position'] = posm 
    varmerged['type_sv'] = indel
    for e in listVar :
        freqtot += e['frequence']
        proftot += e['profondeur']
    varmerged['frequence'] = round(freqtot/len(listVar), 4) #fait la moyenne des fréquences des variants mergés
    varmerged['profondeur'] = int(proftot/len(listVar)) #fait la moyenne des profondeurs des variants mergés
    return varmerged

#fonction qui remplit le dictionnaire de variants mergés, ne retourne rien
def remplitDmerged (DM, Dtout, indel) : #DM le dictionnaire mergé, Dtout un dictionnaire de tous les variants, indel le type de variant
    listpos = sorted(list(Dtout[nom_passage][num_echantillon]))
    i = 0
    #parcourt toutes les positions de l'échantillon
    while i < len(listpos) :
        #cas où on merge avec le suivant :
        if i < len(listpos)-1 and listpos[i] - listpos[i+1] <= 10 : 
            listvar = Dtout[nom_passage][num_echantillon][listpos[i]] + Dtout[nom_passage][num_echantillon][listpos[i+1]]
            pos = int((listpos[i]+listpos[i+1])/2)
            i+=2    #saute le suivant car déjà traité
        #cas simple :
        else :
            listvar = Dtout[nom_passage][num_echantillon][listpos[i]]
            pos = listpos[i]
            i+=1
        remplitDict(DM, mergeVar(listvar, pos, indel))




#parcourt tous les fichiers un par un :
for nom_passage in liste_passages :
    #crée les dictionnaires du passage :
    DtoutINS[nom_passage] = {}
    DtoutDEL[nom_passage] = {}
    Dmerged[nom_passage] = {}

    for num_echantillon in liste_echantillons :
        #crée les dictionnaires de l'échantillon :
        DtoutINS[nom_passage][num_echantillon] = {}
        DtoutDEL[nom_passage][num_echantillon] = {}
        Dmerged[nom_passage][num_echantillon] = {}

        #ouverture du fichier à lire:
        vcf = open(f"{nom_passage}-{num_echantillon}.trimed1000.sv_sniffles.vcf", 'r')

        #parcourt le fichier vcf:
        for l in vcf:   #l chaque ligne du fichier
            
            col = l.split() #récupère les colonnes du vcf

            #ignore les métadonnées et les variants contenant des None:
            if l[0]=='#' or col[1]=='1':
                continue

            #crée un dictionnaire avec les informations du variant
            dvariant = {"passage" : nom_passage, "echantillon": num_echantillon}
            if num_echantillon < 6 :
                dvariant["type_stress"] = "froid"
            else :
                dvariant["type_stress"] = "chaud"
            dvariant["position"] = int(col[1])
            dvariant["seq_reference"] = col[3]
            dvariant["seq_alternative"] = col[4]
            info={}
            for e in col[7].split(';'):
                tmp = e.split('=')
                if tmp[0] == "SVTYPE":
                    info[tmp[0]] = tmp[1]
                if tmp[0] == "AF":
                    info[tmp[0]] = tmp[1]
            dvariant["type_sv"] = info["SVTYPE"]
            dvariant["frequence"] = float(info["AF"])
            dvariant["profondeur"] = int(col[9].split(':')[2]) + int(col[9].split(':')[3])

            #applique les filtres
            if dvariant["frequence"] < seuil_freq or dvariant["profondeur"] < seuil_profondeur :
                continue
            
            #écrit le variant dans le fichier brut
            ecritVariant(dvariant, fichierBrut)

            #remplit le dictionnaire de variants
            if dvariant["type_sv"] == 'INS' :
                remplitDict(DtoutINS, dvariant)
            else:
                remplitDict(DtoutDEL, dvariant)

        #fermeture du fichier à lire:
        vcf.close()

        #crée le dictionnaire Dmerged pour l'échantillon parcourut
        remplitDmerged(Dmerged, DtoutDEL, 'DEL')
        remplitDmerged(Dmerged, DtoutINS, 'INS')
        #et l'écrit dans fichierMerge
        for pos in Dmerged[nom_passage][num_echantillon] :
            for var in Dmerged[nom_passage][num_echantillon][pos] :
                ecritVariant(var, fichierMerge)




#fermeture du fichier de sortie:
fichierBrut.close()
fichierMerge.close()