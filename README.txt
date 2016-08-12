========================================================================================================================
			ClusterGraph
========================================================================================================================
   ***Python 3 est necessaire a l'utilisation de ClusterGraph****
   ***Numpy et scipy sont deux librairies nécessaires aau fonctionnement de la méthode stats pour les graphs***



Tutoriel pour utiliser ClusterGraph:
La liste des options du programme est a la fin de ce fichier. Sinon, vous pouvez utiliser CLusterGraph.py -h

1- Tout d'abord, pour tous les fichiers fasta sur lesquels vous voulez utiliser ClusterGraph, vous devez d'abord utiliser
   prodigal sur ceux-ci pour identifier les genes dans vos sequences.

Exemple:		./prodigal -i fichier_fasta_input.fasta -a fichier_fasta_output.prodigal

2- Ensuite, si vous voulez former un graph avec tous les fichiers sur lesquels vous avez utilise prodigal, vous devez
   modifier le headers de ces fichiers pour savoir qu'elle gene appartient a quel echantillon.

Exemple: J'ai utilise prodigal sur 12a_pne66, 12b_pne66 et 12c_pne66 et les trois fichiers prodigals se situent dans le
         repertoire du projet ClusterGraph/Data/prodigal.

Voici les commandes a faire

- cd ClusterGraph/Data/prodigal

- for i in $(ls); do sed -i s/\>$i\_/g $i; done.

Vos headers dans vos fichiers fasta vont maintenant avoir la forme: 12a_pne66.embl.fasta_x, 12b_pne24.embl.fasta_y,
46_pne66.embl.fasta_z.


3- Ensuite, il suffit de concatene les trois fichiers pour lesquels nous avons modifie les headers pour pouvoir
   ensuite lancer cd-hit sur celui-ci.

Exemple:		 cat *pne* > cat_prodigal.fasta


4- Nous devons maintenant utilise cd-hit pour former nos clusters de genes.

Exemple:		./cd-hit -i cat_prodigal.fasta -o cat_prodigal-cd-hit.fasta -d 0

**L'option -d est important ici pour que cd-hit ne split pas le nom des contigs au premier '_' qu'il rencontre dans
les headers. Nous voulons, par exemple, 12a_pne66_contig-12_1 comme nom de gene et non 12a_1.

5- Lorsque cd-hit termine, il retourne deux fichier: le fichier necessaire parmi ces deux fichiers est le .clstr.
   Ce fichier est le fichier d'entre de ClusterGraph.


6-La demarche pour visualiser votre graph avec Cytoscape ou en html est la suivante:

	6.1- Il faut load le graph avec l'option -i en donnant notre fichier .clstr

	6.2- Il ne faut maintenant qu'utiliser les options -cyto et -j.

	exemple:  python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr -cyto True -j True

Le fichier pour cytoscape est dans le repertoire cytoscape du projet.
Le fichier pour html est dans le repertoire javascript/hcls-dataset-description-master/type-graphs-html du projet.

7-La demarche pour visualise un sous-graph avec des chemins de longueur 10 a partir du cluster 92, par exemple,
  est la suivante:

	7.1- Il faut load le graph avec l'option -i en donnant notre fichier .clstr

		python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr

	7.2- Il faut ensuite utiliser la methode list_of_paths avec l'option -lop qui trouve la liste des chemins
	     du cluster de longueur 10 entre en parametre. On donne le nom du cluster a l'option -lop et la longueur des
	     chemins a l'aide de l'option -x.

		ex: -lop "Cluster 92" -x 10

		***Si vos noms de clusters possede des espaces, il est important de mettre des guillemets autour
		   du nom de votre cluster comme ci-haut.***

	7.3- Ensuite, on cree un sous graph a l'aide de la list_of_paths obtenus grâce aux options -lop et -x avec
	     l'option -sg.
		 ex: -sg True

	7.4- Il ne reste maintenant qu'a utiliser les options -cyto et -j.


	La commande complete est la suivante:


		python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr -lop Cluster 92 -x 10 -sg True -cyto True -j True

Le fichier pour cytoscape est dans le repertoire cytoscape du projet.
Le fichier pour html est dans le repertoire javascript/hcls-dataset-description-master/type-graphs-html du projet.

8- Il existe aussi l'option -s qui est facultative.
	L'option -s ressort quelques statistiques concernant la graph.
		a) Pour tous les noeux du graph, elle ressort le nombre de genes du cluster qui possede le plus de genes.
		b) Pour tous les noeux du graph, elle ressort le nombre de genes du cluster qui possede le moins de genes.
		c) Le nombre de noeux dans le graph.
		d) Le nombre de fois qu'on retrouve 2 fois ou plus deux gene provenenant du même echantillon dans le même cluster.


* Si vous voulez etudier plusieurs genes et que vous possedez un fichier pour lequel la premiere colonne correspond a la
  liste des genes a etudier,vous pouvez passer l'etape 10.

9- etudier plusieurs genes en particulier. Le programme possede un option -xml. Cet option permet de creer des fichiers xml
   lisible dans excel à l'aide de donnée ayant la forme suivante :

Sample_P4J0-Assembly.fa_contig-10000091	27	 16615 	 17580 	 -1 	blaCfxA	beta-lactam	InactivatingEnzyme	82.61	100.00	100.00
Sample_P4J0-Assembly.fa_contig-1251000050	314	 359361 	 360254 	 -1 	blaHGH-1	beta-lactam	InactivatingEnzyme	76.85	98.35	100.00
Sample_P4J0-Assembly.fa_contig-1260000075	99	 109927 	 110886 	 1 	blaHGD-1	beta-lactam	InactivatingEnzyme	99.69	100.00	100.00
Sample_P4J0-Assembly.fa_contig-13000109	6	 3678 	 4523 	 -1 	blaHGH-1	beta-lactam	InactivatingEnzyme	77.30	93.07	100.00
Sample_P4J0-Assembly.fa_contig-2000072	5	 4300 	 5592 	 1 	blaHGH-1	beta-lactam	InactivatingEnzyme	85.81	100.00	70.30
Sample_P4J0-Assembly.fa_contig-26000066	2	 196 	 1191 	 1 	blaCfxA6	beta-lactam	InactivatingEnzyme	98.80	100.00	100.00
Sample_P4J0-Assembly.fa_contig-380000102	3	 2159 	 3133 	 -1 	blaR1	beta-lactam	TranscriptionalRegulator	40.06	54.10	97.54
Sample_P4J0-Assembly.fa_contig-460000086	1	 3 	 1343 	 1 	blaHGH-1	beta-lactam	InactivatingEnzyme	55.29	96.70	65.55
Sample_P4J0-Assembly.fa_contig-562000075	17	 17392 	 18357 	 1 	blaCfxA	beta-lactam	InactivatingEnzyme	82.61	100.00	100.00
Sample_P4J0-Assembly.fa_contig-826000017	1	 1 	 750 	 1 	blaHGH-1	beta-lactam	InactivatingEnzyme	72.18	81.85	99.20
Sample_P4J0-Assembly.fa_contig-89000009	1	 1 	 600 	 1 	blaHGH-1	beta-lactam	InactivatingEnzyme	65.66	65.35	99.00
Sample_P4J7-Assembly.fa_contig-1000054	11	 12655 	 13581 	 -1 	blaHOA-1	beta-lactam	InactivatingEnzyme	40.65	94.56	89.97
Sample_P4J7-Assembly.fa_contig-103	6	 4645 	 5778 	 -1 	omp36	beta-lactam	MutantPorin	70.05	96.81	96.30
Sample_P4J7-Assembly.fa_contig-12000126	2	 2004 	 2828 	 1 	blaCfxA	beta-lactam	InactivatingEnzyme	99.64	85.40	100.00
Sample_P4J7-Assembly.fa_contig-120	15	 13123 	 16272 	 -1 	acrB	acriflavin,aminoglycoside,beta-lactam,glycylcycline,macrolide   Efflux	Efflux	99.90	100.00	100.00
Sample_P4J7-Assembly.fa_contig-120	16	 16295 	 17488 	 -1 	acrA	acriflavin,aminoglycoside,beta-lactam,glycylcycline,macrolide	Efflux	99.50	100.00	100.00

Ces donnees ont ete obtenus a l'aide de la commande suivante:
grep beta /rap/nne-790-ab/projects/Project_CQDM2/Reads-Metagenomics/Sample_P4J7-Assembly/Prodigal/blast/MERGEM-RG.blast.sum.tsv

--->Un fichier possédant la forme suivante est aussi un fichier compatible(la première colonne correspond au gene):

Sample_P4J0-Assembly.fa_contig-10000091_27
Sample_P4J0-Assembly.fa_contig-1251000050_9
Sample_P4J0-Assembly.fa_contig-1260000075_2
Sample_P4J0-Assembly.fa_contig-13000109_6
Sample_P4J0-Assembly.fa_contig-200007_5
Sample_P4J0-Assembly.fa_contig-26000066_2
Sample_P4J0-Assembly.fa_contig-380000102_13
Sample_P4J0-Assembly.fa_contig-460000086_34
Sample_P4J0-Assembly.fa_contig-562000075_5
Sample_P4J0-Assembly.fa_contig-826000017_7
Sample_P4J0-Assembly.fa_contig-89000009_1
Sample_P4J7-Assembly.fa_contig-1000054_11
Sample_P4J7-Assembly.fa_contig-103_6
Sample_P4J7-Assembly.fa_contig-12000126?2
Sample_P4J7-Assembly.fa_contig-120_15
Sample_P4J7-Assembly.fa_contig-120_16


Pour utiliser l'option -xml, il suffit donc de stocker les donnees obtenues a l'aide du grep dans un fichier.
Exemple:
grep beta /rap/nne-790-ab/projects/Project_CQDM2/Reads-Metagenomics/Sample_P4J7-Assembly/Prodigal/blast/MERGEM-RG.blast.sum.tsv > ClusterGraph/Data/XML/Beta_Lactam_Fred_P4J0-7-90

Voici un exemple d'utilisation de l'option -xml. Nous allons donc donner le fichier Beta_Lactam_Fred_P4J0-7-90 au programme
de cette façon en ajoutant l'option -xml. Cet option permet de creer un fichier qui peut s'ouvrir dans excel qui compare
les chemins des differents genes qui etaient presents dans le fichier Beta_Lactam_Fred_P4J0-7-90:

	python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr  -xml ClusterGraph/Data/XML/Beta_Lactam_Fred_P4J0-7-90 -x 10

Pour spécifier la longueur des chemins pour le fichier xml, enmcore une fois utiliser l'option -x.

	Si vous possedez un fichier pour lequel la premiere colonne correspond a la liste des genes a etudier, remplacez Beta_Lactam_Fred_P4J0-7-90 par votre fichier avec son chemin.
	De plus, dans ce cas, vous pouvez utiliser l'option -mn qui vous permet de creer plusieurs fichier cytoscape pour tous les genes qui se trouvent dans votre fichier.
	Ces fichiers cytoscape se retrouvent dans le repertoire many_genes_cytoscape dans le repertoire Cytoscape.

	python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr  -xml votre_fichier -x 10

	cytoscape:
	python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr -mn fichier_contenant_vos_genes -x 10



10- Il est aussi possible de sauvegarder le graph apres l'avoir load, puis de le reloader par la suite a partir du fichier de sauvegarde.

Sauvegarder:
	python3 ClusterGraph.py -save sauvegarde_graph.txt

Reload:
	python3 ClusterGraph.py -r sauvegarde_graph.txt



---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
			ClusterGraph.py -h
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  -h, --help  show this help message and exit

  -i I        Input file.Le fichier en input doit correspondre au fichier
              .clstr de cd-hit.

  -find FIND  L'input doit être le nom d'un gene qui est present dans le
              graph. La fonction va alors retourner le cluster respectif du
              gene.

  -lop LOP    La fonction retourne tous les chemins de longueur "x" partant du
              nom du cluster que vous devez donner en parametetre. Utilisez -x
              pour donner la longueur des chemins desires.

  -x X        Longueur des chemin desire lors de l'utilisation de -lop ou -mn
              ou les deux.

  -sg SG      Construit un sous-graph a partir des donnees reçues par -lop et
              -x| True ou False.

  -cyto CYTO  Permet de visualiser le graph dans cytoscape a l'aide du fichier
              cytoscape.txtCelui-ci doit se trouve dans le repertoire
              Cytoscape du projet.| True ou False.

  -j J        Permet de visualiser le graph dans une page html de façon
              interactive.Il ne faut qu'ouvrir le fichier index.html dans
              javascript/hcls-dataset-description-master/type-graphs-html.|
              True ou False.

  -xml XML    Creer un fichier xml qui peut être ouvert dans Microsoft Excel.|
              True ou False.

  -s S        Affiche differentes statistiques sur le graph.

  -save SAVE  Sauvegarde le graph dans le fichier donne en parametre.

  -r R        Load le graph a partir du fichier de sauvegarde deonne en
              parametre

  -mn MN      Creer des fichiers cytoscape pour plusieurs genes a la
              fois.Utilisez -x quand vous utilisez -mn pour specifier la
              longueur des chemins.











