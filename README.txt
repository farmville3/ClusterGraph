========================================================================================================================
			ClusterGraph
========================================================================================================================
   ***Python 3 est necessaire à l'utilisation de ClusterGraph****



Tutoriel pour utiliser ClusterGraph:
La liste des options du programme est à la fin de ce fichier. Sinon, vous pouvez utiliser CLusterGraph.py -h

1- Tout d'abord, pour tous les fichiers fasta sur lesquels vous voulez utiliser ClusterGraph, vous devez d'abord utiliser
   prodigal sur ceux-ci pour identifier les genes dans vos sequences.

Exemple:		./prodigal -i fichier_fasta_input.fasta -a fichier_fasta_output.prodigal 

2- Ensuite, si vous voulez former un graph avec tous les fichiers sur lesquels vous avez utilise prodigal, vous devez
   modifier le headers de ces fichiers pour savoir qu'elle gene appartient à quel echantillon.

Exemple: J'ai utilise prodigal sur 12a_pne66, 12b_pne66 et 12c_pne66 et les trois fichiers prodigals se situent dans le
         repertoire du projet ClusterGraph/Data/prodigal.

Voici les commandes à faire

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

7-La demarche pour visualise un sous-graph avec des chemins de longueur 10 à partir du cluster 92, par exemple,
  est la suivante:

	7.1- Il faut load le graph avec l'option -i en donnant notre fichier .clstr

		python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr

	7.2- Il faut ensuite utiliser la methode list_of_paths avec l'option -lop qui trouve la liste des chemins
	     du cluster de longueur 10 entre en parametre. On donne le nom du cluster à l'option -lop et la longueur des
	     chemins à l'aide de l'option -x.

		ex: -lop "Cluster 92" -x 10

		***Si vos noms de clusters possede des espaces, il est important de mettre des guillemets autour
		   du nom de votre cluster comme ci-haut.***

	7.3- Ensuite, on cree un sous graph à l'aide de la list_of_paths obtenus grâce aux options -lop et -x avec
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


* Si vous voulez etudier plusieurs genes et que vous possedez un fichier pour lequel la premiere colonne correspond à la
  liste des genes à etudier,vous pouvez passer l'etape 10.

9- etudier plusieurs genes en particulier. Le programme possede un option -g. Cet option permet de convertir les donnees
   possedant cette forme :

contig-1000054	11	 12655 	 13581 	 -1 	blaHOA-1	beta-lactam	InactivatingEnzyme	40.65	94.56	89.97
contig-103	6	 4645 	 5778 	 -1 	omp36	beta-lactam	MutantPorin	70.05	96.81	96.30
contig-12000126	2	 2004 	 2828 	 1 	blaCfxA	beta-lactam	InactivatingEnzyme	99.64	85.40	100.00
contig-120	15	 13123 	 16272 	 -1 	acrB	acriflavin,aminoglycoside,beta-lactam,glycylcycline,macrolide   Efflux	Efflux	99.90	100.00	100.00
contig-120	16	 16295 	 17488 	 -1 	acrA	acriflavin,aminoglycoside,beta-lactam,glycylcycline,macrolide	Efflux	99.50	100.00	100.00
contig-129000083	1	 1 	 702 	 1 	blaHGH-1	beta-lactam	InactivatingEnzyme	62.88	75.58	97.86
contig-149000032	2	 380 	 1375 	 1 	blaCfxA6	beta-lactam	InactivatingEnzyme	98.80	100.00	100.00


Ces donnees ont ete obtenus à l'aide de la commande suivante:
grep beta /rap/nne-790-ab/projects/Project_CQDM2/Reads-Metagenomics/Sample_P4J7-Assembly/Prodigal/blast/MERGEM-RG.blast.sum.tsv

Pour utiliser l'option -g, il suffit donc de stocker les donnees obtenues à l'aide du grep dans un fichier.
Exemple:
grep beta /rap/nne-790-ab/projects/Project_CQDM2/Reads-Metagenomics/Sample_P4J7-Assembly/Prodigal/blast/MERGEM-RG.blast.sum.tsv > ClusterGraph/Data/XML/Beta_Lactam_Fred_P4J0-7-90

10- Nous allons donc donner le fichier test.txt au programme de cette façon en ajoutant l'option -xml. Cet option
    creer un fichier qui peut s'ouvrir dans excel qui compare les chemins des differents genes qui etaient presents dans le fichier test.txt:
	python3 ClusterGraph.py -g ClusterGraph/Data/XML/Beta_Lactam_Fred_P4J0-7-90 -xml True

	Si vous possedez un fichier pour lequel la premiere colonne correspond à la liste des genes a etudier, remplacez test.txt par votre fichier avec son chemin.
	De plus, dans ce cas, vous pouvez utiliser l'option -mn qui vous permet de créer plusieurs fichier cytoscape pour tous les gènes qui se trouvent dans votre fichier.

	python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr  -g votre_fichier -xml True -y 10

	cytoscape:
	python3 ClusterGraph.py -i cat_prodigal-cd-hit.fasta.clstr -mn fichier_contenant_vos_genes -x 10



11- Il est aussi possible de sauvegarder le graph apres l'avoir load, puis de le reloader par la suite à partir du fichier de sauvegarde.

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

  -find FIND  L'input doit être le nom d'un gène qui est présent dans le
              graph. La fonction va alors retourner le cluster respectif du
              gène.

  -lop LOP    La fonction retourne tous les chemins de longueur "x" partant du
              nom du cluster que vous devez donner en parametètre. Utilisez -x
              pour donner la longueur des chemins désirés.

  -x X        Longueur des chemin désiré lors de l'utilisation de -lop ou -mn
              ou les deux.

  -sg SG      Construit un sous-graph à partir des données reçues par -lop et
              -x| True ou False.

  -cyto CYTO  Permet de visualiser le graph dans cytoscape à l'aide du fichier
              cytoscape.txtCelui-ci doit se trouvé dans le répertoire
              Cytoscape du projet.| True ou False.

  -j J        Permet de visualiser le graph dans une page html de façon
              interactive.Il ne faut qu'ouvrir le fichier index.html dans
              javascript/hcls-dataset-description-master/type-graphs-html.|
              True ou False.

  -g G        Convertit les fichiers obtenus à l'aide de greps en ligne de
              commande de façonà ce que les données soient convertible dans un
              fichier .xml

  -xml XML    Créer un fichier xml qui peut être ouvert dans Microsoft Excel.|
              True ou False.

  -y Y        Longueur des chemin désiré lors de l'utilisation de -xml.

  -s S        Affiche différentes statistiques sur le graph.

  -save SAVE  Sauvegarde le graph dans le fichier donné en paramètre.

  -r R        Load le graph à partir du fichier de sauvegarde deonné en
              paramètre

  -mn MN      Créer des fichiers cytoscape pour plusieurs gènes à la
              fois.Utilisez -x quand vous utilisez -mn pour spécifier la
              longueur des chemins.











