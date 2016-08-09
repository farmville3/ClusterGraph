==================================================================================================
			                      ClusterGraph
==================================================================================================

Tutoriel pour utiliser ClusterGraph:
La liste des options du programme est à la fin de ce fichier. Sinon, vous pouvez utiliser CLusterGraph.py -h

1- Tout d'abord, pour tous les fichiers fasta sur lesquels vous voulez utiliser ClusterGraph, vous devez d'abord utiliser prodigal sur ceux-ci pour identifier les gènes dans vos séquences.

Exemple:		./prodigal -i fichier_fasta_input.fasta -a fichier_fasta_output.prodigal 

2- Ensuite, si vous voulez former un graph avec tous les fichiers sur lesquels vous avez utilisé prodigal, vous devez modifier le headers de ces fichiers pour savoir qu'elle gène appartient à quel échantillon.

Exemple: J'ai utilisé prodigal sur 12a_pne66, 12b_pne66 et 12c_pne66 et les trois fichiers prodigals se situent dans le répertoire ~/prodigal_data.

Voici les commandes à faire

- cd ~/prodigal_data

- for i in $(ls); do sed -i s/\>$i\_/g $i; done.

Vos headers dans vos fichiers fasta vont maintenant avoir la forme: 12a_pne66_contig-x, 12b_pne66_contig-y,12c_pne66_contig-z.


3- Ensuite, il suffit de concatené les trois fichiers pour lesquels nous avons modifié les headers pour pouvoir ensuite lancer cd-hit sur celui-ci.

Exemple:		 cat *pne* > all__pne66.prodigal


4- Nous devons maintenant utilisé cd-hit pour former nos clusters de gènes.

Exemple:		./cd-hit -i all_pne66.prodigal -o all_pne66_cd-hit.prodigal -d 0

**L'option -d est important ici pour que cd-hit ne split pas le nom des contigs au premier '_' qu'il rencontre dans les headers.
Nous voulons, par exemple, 12a_pne66_contig-12_1 comme nom de gène et non 12a_1.

5- Lorsque cd-hit termine, il retourne deux fichier: le fichier nécessaire parmi ces deux fichiers est le .clstr. Ce fichier est le fichier d'entré de ClusterGraph.


6-La démarche pour visualiser votre graph avec Cytoscape ou en html est la suivante:

	6.1- Il faut load le graph avec l'option -i en donnant notre fichier .clstr

	6.2- Il ne faut simplement qu'utiliser les options -cyto et -j.

	exemple:  python3 ClusterGraph.py -i all_pne66_cd-hit.prodigal.clstr -cyto True -j True

Le fichier pour cytoscape est dans le répertoire cytoscape du projet.
Le fichier pour html est dans le répertoire javascript/hcls-dataset-description-master/type-graphs-html du projet.

7-La démarche pour visualisé un sous-graph avec des chemins de longueur 10 à partir du cluster 92, par exemple, est la suivante:

	7.1- Il faut load le graph avec l'option -i en donnant notre fichier .clstr

		python3 ClusterGraph.py -i all_pne66_cd-hit.prodigal.clstr

	7.2- Il faut ensuite utiliser la méthode list_of_paths avec l'option -lop qui trouve la liste des chemins du cluster de longueur 10 entré en paramètre.
		On donne le nom du cluster à l'option -lop et la longueur des chemins à l'aide de l'option -x.

		ex: -lop "Cluster 92" -x 10

		***Si vos noms de clusters possède des espaces, il est important de mettre des guillemets autour du nom de votre cluster comme ci-haut.***

	7.3- Ensuite, on créer un sous graph à l'aide de la list_of_paths obtenus grâce aux options -lop et -x avec l'option -sg.
		ex: -sg True

	7.4- Il ne reste maintenant qu'a utiliser les options -cyto et -j.

	
	La commande complète est la suivante:

	
		python3 ClusterGraph.py -i all_pne66_cd-hit.prodigal.clstr -lop Cluster 92 -x 10 -sg True -cyto True -j True

Le fichier pour cytoscape est dans le répertoire cytoscape du projet.
Le fichier pour html est dans le répertoire javascript/hcls-dataset-description-master/type-graphs-html du projet.

8- Il existe aussi l'option -s qui est facultative.
	L'option -s ressort quelques statistiques concernant la graph.
		a) Pour tous les noeux du graph, elle ressort le nombre de gènes du cluster qui possède le plus de gènes.
		b) Pour tous les noeux du graph, elle ressort le nombre de gènes du cluster qui possède le moins de gènes.
		c) Le nombre de noeux dans le graph.
		d) Le nombre de fois qu'on retrouve 2 fois ou plus deux gène provenenant du même échantillon dans le même cluster.

9- Étudier plusieurs gènes en particulier. Le programme possède un option -g. Cet option permet de convertir les données possédant cette forme :

contig-1000054	11	 12655 	 13581 	 -1 	blaHOA-1	beta-lactam	InactivatingEnzyme	40.65	94.56	89.97
contig-103	6	 4645 	 5778 	 -1 	omp36	beta-lactam	MutantPorin	70.05	96.81	96.30
contig-12000126	2	 2004 	 2828 	 1 	blaCfxA	beta-lactam	InactivatingEnzyme	99.64	85.40	100.00
contig-120	15	 13123 	 16272 	 -1 	acrB	acriflavin,aminoglycoside,beta-lactam,glycylcycline,macrolide   Efflux	Efflux	99.90	100.00	100.00
contig-120	16	 16295 	 17488 	 -1 	acrA	acriflavin,aminoglycoside,beta-lactam,glycylcycline,macrolide	Efflux	99.50	100.00	100.00
contig-129000083	1	 1 	 702 	 1 	blaHGH-1	beta-lactam	InactivatingEnzyme	62.88	75.58	97.86
contig-149000032	2	 380 	 1375 	 1 	blaCfxA6	beta-lactam	InactivatingEnzyme	98.80	100.00	100.00


Ces données ont été obtenus à l'aide de la commande suivante: grep beta /rap/nne-790-ab/projects/Project_CQDM2/Reads-Metagenomics/Sample_P4J7-Assembly/Prodigal/blast/MERGEM-RG.blast.sum.tsv

Pour utiliser l'option -g, il suffit donc de stocker les données obtenues à l'aide du grep dans un fichier.
Exemple:	grep beta /rap/nne-790-ab/projects/Project_CQDM2/Reads-Metagenomics/Sample_P4J7-Assembly/Prodigal/blast/MERGEM-RG.blast.sum.tsv > test.txt

Nous allons donc donner le fichier test.txt au programme de cette façon en ajoutant l'option -xml. Cet option créer un fichier qui peut s'ouvrir dans excel qui compare les chemins des différents gènes qui étaient présents dans le fichier test.txt:
	python3 ClusterGraphy.py -g test.txt -xml True


10- Il est aussi possible de sauvegarder le graph après l'avoir load, puis de le reloader par la suite à partir du fichier de sauvegarde.

Sauvegarder:
	python3 ClusterGraph.py -save sauvegarde_graph.txt

Reload:
	python3 ClusterGraph.py -r sauvegarde_graph.txt



------------------------------------------------------------------------------------------------------
			            ClusterGraph.py -h
------------------------------------------------------------------------------------------------------
optional arguments:
  -h, --help  show this help message and exit

  -i I        Input file.Le fichier en input doit correspondre au fichier
              .clstr de cd-hit.

  -find FIND  L'input doit être le nom d'un gène qui est présent dans le
              graph. La fonction va alors retourner le cluster respectif du
              gène.

  -lop LOP    La fonction retourne tous les chemins de longueur "x" partant du
              nom du cluster que vous devez donner en parametètre. Utilisez -x
              pour donner la longueur des chemins désirés.
  -x X        Longueur des chemin désiré lors de l'utilisation de -lop.
  -sg SG      Construit un sous-graph à partir des données reçues par -lop et
              -x

  -cyto CYTO  Permet de visualiser le graph dans cytoscape à l'aide du fichier
              cytoscape.txt Celui-ci doit se trouvé dans le répertoire
              Cytoscape du projet.

  -j J        Permet de visualiser le graph dans une page html de façon
              interactive.Il ne faut qu'ouvrir le fichier index.html dans
              javascript/hcls-dataset-description-master/type-graphs-html
  -g G        Convertit les fichiers obtenus à l'aide de greps en ligne de
              commande de façonà ce que les données soient convertible dans un
              fichier .xml

  -xml XML    Créer un fichier xml qui peut être ouvert dans Microsoft Excel.
  -y Y        Longueur des chemin désiré lors de l'utilisation de -xml.
  -s S        Affiche différentes statistiques sur le graph.

  -save SAVE  Sauvegarde le graph dans le fichier donné en paramètre.
  -r R        Load le graph à partir du fichier de sauvegarde deonné en
              paramètre









