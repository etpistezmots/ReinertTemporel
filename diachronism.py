#!/usr/bin/env python
# -*- coding: utf-8 -*-

# fait par Max Beligné

from math import sqrt

# Vesion Python du code suivant
# https://github.com/nicolasdugue/istex/blob/master/JAVA/CSR/src/model/diachronism/LabelDiachronism.java
# format des fichiers d'entrée :
# https://github.com/nicolasdugue/istex/tree/master/JAVA/Expe_Scientometrics
# These files are results from labeling of the clustering results of two distinct periods

# 2 Choix de jeux test dispo : "Nicolas" ou "Max"
ChoixJeuTest = "Max"

if ChoixJeuTest == "Nicolas":
    periode1 = __file__[:-15] + "/Data/00-0350gnge-tn-DD-ENL.fmgs"
    periode2 = __file__[:-15] + "/Data/96-99-F43gnge-tn-DD-ENL.fmgs"

if ChoixJeuTest == "Max":
    periode1 = __file__[:-15] + "/Data/spe1newpourcent.fmgs"
    periode2 = __file__[:-15] + "/Data//spe2newpourcent.fmgs"


class diachronism:

    def __init__(self, periode1, periode2, env="JusteActif"):

        ##########  VARIABLES #################

        # les dicos me permettent de faire des appels communs
        # et de marquer dans les prints de quoi il s'agit grace à l'argumen "name"
        self.P = {"P1": {}, "P2": {}}
        self.P1 = {"data" : self.P["P1"], "name" : "P1"}
        self.P2 = {"data" : self.P["P2"], "name" : "P2"}

        self.PTsachantS = {"P1VersP2": [], "P2VersP1": []}
        # je ne suis pas convaincu par la notation suivante
        # qui mélange les périodes et les termes "source" et "cible"
        # en fait, on va toujours de la source vers la cible
        # mais bon en ref à l'article http://lodel.irevues.inist.fr/isj/?id=390:
        self.Ps_t = self.PTsachantS["P1VersP2"]
        self.Pt_s = self.PTsachantS["P2VersP1"]
        # il faut dans le cas de cette notation comprendre source = P1 et cible = P2

        self.PA = {"P1" : [], "P2" : []}
        # idem precédent
        self.PA_s = self.PA["P1"]
        self.PA_t = self.PA["P2"]

        self.A = {"P1": 0, "P2": 0}
        # idem
        self.A_s = self.A["P1"]
        self.A_t = self.A["P1"]

        self.D = {"P1": 0, "P2": 0}
        # id
        self.D_s = self.D["P1"]
        self.D_t = self.A["P1"]

        self.ListMatch = []


        ################# OBTENTION DES DONNEES ###################
        #  variable : self.P

        self.GetLabelAndScore(periode1, "P1")
        self.GetLabelAndScore(periode2, "P2")


        ############ MATCHING (Cf http://lodel.irevues.inist.fr/isj/?id=390) ##########

            # Compute P(T|S) et P(S|T)  --> Equation 13
            # variable : self.PTsachantS


        self.ComputeProbTKnowingS(self.P1, self.P2, "P1VersP2")
        self.ComputeProbTKnowingS(self.P2, self.P1, "P2VersP1")

            # Compute PA(S) et PA(T) --> Equation 14
            # variable : self.PA

        self.ComputeLocalAverage(self.P1, self.PTsachantS["P1VersP2"], env)
        self.ComputeLocalAverage(self.P2, self.PTsachantS["P2VersP1"], env)

            # Compute Global Average A(S) et A(T) --> Equation 15
            # variable : self.A

        self.ComputeGlobalAverage(self.P1)
        self.ComputeGlobalAverage(self.P2)

            # Compute standard Deviation
            # variable : self.D

        self.ComputeStandardDeviation(self.P1)
        self.ComputeStandardDeviation(self.P2)

            # Match cluster
            # variable : self.ListMatch
        self.MatchCluster(self.P1, self.P2, self.PTsachantS["P1VersP2"], self.PTsachantS["P2VersP1"], self.PA["P1"], self.PA["P2"], self.A["P1"], self.A["P2"], self.D["P1"], self.D["P2"],)



    ###############  FONCTIONS DE L'INIT ###################

    def GetLabelAndScore(self, pathdata, indic):

        '''
        En sortie, renvoie un dico de dico sous la forme suivante
        {"P1" : {"label1" : score1, "label2" : score2,...}, "nom ensemble 2" :....  }
        qui sera lui même stocké dans le dico self.P
        '''

        NomPartie = ""

        # lit le fichier où sont stocké les données
        with open(pathdata) as f:
            content = f.readlines()

        for line in content:
            # Si c'est la marque d'un nouvel ensemble
            # Ex "G0-0	25	19" avec des tabulations comme séparateur
            if line[0] == "G":
                # recup son nom"
                EnsInfo = line.split("\t")
                NomPartie = EnsInfo[0][1:].split("-")[0] + "_" + indic

                # Crée une nouvel entrée dans le dico
                self.P[indic][NomPartie] = {}

            # Ex "	industrie	5488" mais aussi espace entre les ensembles (d'où le "and line[0:1]!="\n")
            if line[0] != "G" and line[0:1]!="\n":
                mot, score = line.split("\t")[1], line.split("\t")[2][:-1]

                # si correspond à un nouveau mot, insertion dans le dico avec le score correspondant
                if mot != "":
                    self.P[indic][NomPartie][mot] = float(score)



    def ComputeProbTKnowingS(self, AllEnsSource, AllEnsCible, indic):
        '''
        Indic : permet de clarifier et faire appel à dico pour stocker résultat

        En sortie, liste de dico {"cible", "source", "proba de l'elt cible sachant l'elt source")
        Seuls les elts où il y a eu une intersection entre elt source et elt cible sont marqués
        Pour les autres, la proba est égal à O car il n'y a rien au numérateur
        '''


        # boucle sur l'ensemble combinaision ensemble source- ensemble cible
        # Sachant S, compute P(t|s)

        print("\n")
        print("--- DETAIL COMPUTE PROBABILITE DE " + AllEnsCible["name"] + " SACHANT "  +AllEnsSource["name"] + "-------------")

        # sachant s
        for eltS in AllEnsSource["data"]:

            # compute t
            for eltT in AllEnsCible["data"]:

                print(eltS, eltT)

                # recherche étiquette intersection element source - element cible
                EltIntS = AllEnsSource["data"][eltS].keys()
                EltIntT = AllEnsCible["data"][eltT].keys()
                intersection = EltIntS & EltIntT
                print(intersection)

                if len(intersection) > 0:

                    # Calcul du numerateur = somme score des étiquettes instersection (sachant S) dans la destination
                    numerateur = 0
                    for eltI in intersection:
                        numerateur += AllEnsCible["data"][eltT][eltI]
                    print(numerateur)

                    # Calcul de dénominateur = somme score étiquette dans destination

                    denominateur = 0
                    for eltLabelT in AllEnsCible["data"][eltT]:
                        denominateur += AllEnsCible["data"][eltT][eltLabelT]
                    print(denominateur)

                    # Calul de P(t|s)
                    if denominateur != 0:
                        ProbTS = numerateur / denominateur
                    else:
                        ProbTS = 0

                    print("La ProbTS de " + eltT + " sachant " + eltS + " est " + str(ProbTS))
                    self.PTsachantS[indic].append({"cible" : eltT, "source": eltS, "proba" : ProbTS})

        print("-- FIN DETAIL COMPUTE PROBABILITE DE " + AllEnsCible["name"] + " SACHANT "  +AllEnsSource["name"] + "-------------")

        print("\n")
        print("-- RESULTAT COMPUTE PROBABILITE DE " + AllEnsCible["name"] + " SACHANT "  +AllEnsSource["name"] + "-------------")

        print(self.PTsachantS[indic])



    def ComputeLocalAverage(self, Ens, ProbaResultEns, testenv):

        '''Variable testenv car question sur quoi on divise
        que les clusters actifs : implémentation actuelle de Nicolas. Variable : JusteActif
        Quand beaucoup de clusters comme dans le jeu Test Nicolas Préférable
        
        Quand peu de clusters comme dans mon jeu de test perso
        Peut se justifier : voir cas suivant
        
        - Le cluster 0 en période 1 active seulement le cluster 5 en période 2 à hauteur de 0.31
        - Lecluster 5 en période 2 active réciproquement le cluster 0 en période 1 à hauteur de 0.43
        mais aussi un tout petit peu (pour un seul mot) le cluster 3 en période 1 à hauteur de 0.04. Du coup, ça fait chuter sa moyenne à 0,235.
        - Comme on prend la moyenne interne, ce tout petit mot à un impact énorme'''


        print("\n")
        print("--- DETAIL COMPUTE LOCAL AVERAGE " + Ens["name"] + "  -------------")

        for elt in Ens["data"]:
            print(elt)
            # Quand dans les résultats précédents l'element source est égal à elt
            ListeActivateByElt = [result for result in ProbaResultEns if result["source"] == elt]
            # Dans ce cas là, elt activé par la source, on somme les proba des elts activés
            SommePTsachantS = sum([result["proba"] for result in ListeActivateByElt])

            if testenv == "JusteActif":
                print(SommePTsachantS)
                Env = len(ListeActivateByElt)
                print(Env)
                if Env != 0:
                    ProbaActiviteS = SommePTsachantS / Env
                else:
                    ProbaActiviteS = 0
                print(ProbaActiviteS)
                print("\n")

            # ESSAI ATTENTION si n'est pas égal à valeur par défaut ("JusteActif") :
            # fait automatiquement le calcul suivant !
            else :
                print(SommePTsachantS)
                print(len(Ens["data"]))
                ProbaActiviteS = SommePTsachantS / len(Ens["data"])
                print(ProbaActiviteS)
                print("\n")

            self.PA[Ens["name"]].append({"nom" : elt, "activite" : ProbaActiviteS})

        print("---------- FIN DETAIL COMPUTE LOCAL AVERAGE " + Ens["name"] + "-------------")

        print("\n")
        print("---------- RESULTAT COMPUTE LOCAL AVERAGE " + Ens["name"] + "-----------")
        print("\n")
        print(self.PA[Ens["name"]])


    def ComputeGlobalAverage(self, Ens):

        # recup du resultat précédent des activation provoqués par l'ensemble d'une période
        Activ = sum([elt["activite"] for elt in self.PA[Ens["name"]]])
        NbreEns = len(self.PA[Ens["name"]])
        self.A[Ens["name"]] = Activ / NbreEns
        print("\n")
        print("---------- RESULTAT COMPUTE GLOBAL AVERAGE " + Ens["name"] + "-----------")
        print("\n")
        print(self.A[Ens["name"]])


    def ComputeStandardDeviation(self, Ens):

        observationsP = [elt["activite"] for elt in self.PA[Ens["name"]]]
        MoyenneP = self.A[Ens["name"]]

        VarianceElt = [(elt - MoyenneP) ** 2 for elt in observationsP]
        Variance = sum(VarianceElt) / len(VarianceElt)
        self.D[Ens["name"]] = sqrt(Variance)
        print("\n")
        print("---------- RESULTAT COMPUTE STANDARD DEVIATION SOURCE -----------")
        print("\n")
        print(self.D[Ens["name"]])



    def MatchCluster(self, AllEnsSource, AllEnsCible, AllProba2sachant1, AllProba1sachant2, PA1, PA2, A1, A2, D1, D2):

        # itération sur toutes les combinaisons
        for eltS in AllEnsSource["data"]:
            for eltT in AllEnsCible["data"]:

                # recupère les proba calulées
                Proba2sachant1List = [elt["proba"] for elt in AllProba2sachant1 if (elt["cible"] == eltT and elt["source"]== eltS)]
                Proba1sachant2List = [elt["proba"] for elt in AllProba1sachant2 if (elt["cible"] == eltS and elt["source"] == eltT)]

                # Gère le cas où il n'y a pas de proba (=0)
                if len(Proba2sachant1List)!=0:
                    Proba2sachant1= Proba2sachant1List[0]
                else:
                    Proba2sachant1 = 0

                if len(Proba1sachant2List)!=0:
                    Proba1sachant2 = Proba1sachant2List[0]
                else:
                    Proba1sachant2 = 0


                # récupere le PAs et Pat
                SpecificPA1 = [elt["activite"] for elt in PA1 if elt["nom"]==eltS][0]
                SpecificPA2 = [elt["activite"] for elt in PA2 if elt["nom"] == eltT][0]

                # ensemble de print pour vérification
                if len(Proba2sachant1List) != 0 and len(Proba1sachant2List)!=0:
                    print("Part1 : " + eltS)
                    print("Part2 : " + eltT)
                    print("PPart2SachantPart1 : " + str(Proba2sachant1))
                    print("PPart1SachantPart2 : " + str(Proba1sachant2))
                    print("PA_part1 : " + str(SpecificPA1))
                    print("PA_part2 : " + str(SpecificPA2))
                    print("A_part1 + D_part1 : " + str(A1 + D1))
                    print("A_part2 + D_part2 : " + str(A2 + D2))
                    print("\n")

                # Def du matching : voir http://lodel.irevues.inist.fr/isj/?id=390
                # après l'aquation 15
                if Proba2sachant1 >= SpecificPA1 and (Proba2sachant1 >= A1 + D1) and Proba1sachant2 >= SpecificPA2 and (Proba1sachant2 >= A2 + D2):
                    self.ListMatch.append((eltS,eltT))



############################ TEST APPEL #########################


# EXPE 1
print("\n")
print("-----------EXPE1 DIVISION SUR JUSTE CLUSTER ACTIF--------------")
TestMax = diachronism(periode1,periode2)
print("\n")
print("-----------FIN EXPE 1 DIVISION SUR JUSTE CLUSTER ACTIF-------------")

# EXPE 2
print("\n")
print("-----------EXPE2 DIVISION SUR TOUS LES CLUSTERS--------------")
TestMax2 = diachronism(periode1,periode2,env="Test")
print("\n")
print("-----------FIN EXPE2 DIVISION SUR TOUS LES CLUSTERS--------------")

# PRINT RESULTAT EXPE 1 ET EXPE 2 !
print("\n")
print("RESULT FIN AVEC DIVISION SUR JUSTE CLUSTERS ACTIFS")
print(TestMax.ListMatch)


print("\n")
print("RESULT FIN AVEC DIVISION SUR TOUS LES CLUSTERS")
print(TestMax2.ListMatch)


# Attention la fin n'est pas implémentée (à partir l 270)
# rajout règle pas vraiment justifiée à mon avis
# si 1) = p_t_knowing_s et 2) = p_s_knowing_t
# 1er cas = [1)*0,666 ≤  2)  et    2)*0,666 >=1) ou [ 1)<0,34 et 2) <0,34] ou [ 1)>0,32 et 2)>0,32]
# problème ne couvre pas mathématiquement l’ensemble des cas : fait plus penser à bricolage

