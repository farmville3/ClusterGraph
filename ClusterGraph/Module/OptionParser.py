import argparse
import sys

#----------------------OptionParser----------------------------------------------------------------------------------------------------------------------------------------
class OptionParser():
    def __init__(self, args):
        self.parser = argparse.ArgumentParser(description="Auto phage extractor<eta.")


        ##########         parser       #########################
        self.parser.add_argument('-i', type=str, help='Input file.\nLe fichier en input doit correspondre au fichier .clstr de cd-hit.',default='',required=True)

        self.parser.add_argument('-find', type=str, help="L'input doit être le nom d'un gène qui est présent dans le graph.\n La fonction va alors retourner le cluster respectif du gène.", default='',
                                         required=False)

        self.parser.add_argument('-lop', type=str, help='La fonction retourne tous les chemins de longueur "x" partant du nom du cluster que vous devez donner en parametètre.'
                                                        '\n Utilisez -x pour donner la longueur des chemins désirés.', default='',required=False)

        self.parser.add_argument('-x', type=int, help="Longueur des chemin désiré lors de l'utilisation de -lop.", default=-1,
                                 required=False)

        self.parser.add_argument('-sg', type=bool, help="Construit un sous-graph à partir des données reçues par -lop et -x", default=False,
                                 required=False)

        self.parser.add_argument('-cyto', type=bool, help="Permet de visualiser le graph dans cytoscape à l'aide du fichier cytoscape.txt\n"
                                 'Celui-ci doit se trouvé dans le répertoire Cytoscape du projet. ', default=False,required=False)

        self.parser.add_argument('-j', type=bool, help="Permet de visualiser le graph dans une page html de façon interactive."
                                 "Il ne faut qu'ouvrir le fichier index.html dans javascript/hcls-dataset-description-master/type-graphs-html", default=False,required=False)

        self.parser.add_argument('-g', type=str, help="Convertit les fichiers obtenus à l'aide de greps en ligne de commande de façon\n"
                                                      "à ce que les données soient convertible dans un fichier .xml", default='',
                                 required=False)

        self.parser.add_argument('-xml', type=bool, help="Créer un fichier xml qui peut être ouvert dans Microsoft Excel.", default='',
                                 required=False)

        self.parser.add_argument('-y', type=int, help="Longueur des chemin désiré lors de l'utilisation de -xml.", default=-1,
                                 required=False)

        self.parser.add_argument('-s', type=bool, help="Affiche différentes statistiques sur le graph.", default=True,
                                 required=False)

        self.parser.add_argument('-t', type=bool, help="Affiche le temps écoulé entre le début et la fin de l'appel de la fonction.", default='',
                                 required=False)



        #parse the args...
        self.Arguments = vars(self.parser.parse_args(args))

    #Needed to get the args
    def getArguments(self):
        return self.Arguments

    def getParser(self):
        return self.parser

if __name__ == '__main__':
    print('========================================================================================================')
    print("\t"+"\t"+"\t"+"\t"+"\t"+"\t"+"\t"+"\t"+"\t"+"\t"+"\t"+"\t"+'\033[4m' +'start'+ '\033[0m')
    print('========================================================================================================')

    parser = OptionParser(sys.argv[1:])
    arg = parser.getArguments()
    load_file= arg['i']
    find = arg['find']
    list_of_paths=arg['lop']
    lenght_lop=arg['x']
    sous_graph=arg['sg']
    cytoscape=arg['cyto']
    javacript=arg['j']
    convert_greps_fred=arg['g']
    xml = arg['xml']
    lenght_xml=arg['y']
    stats=arg['s']
    time=arg['t']

    #Load graph
    graph = Graph()

    if load_file!='':
        #graph.load_graph(('/home/saiant01/cat_Sample_P4Jx-Assembly_100.fa.clstr'))

        #P4
        graph.load_graph(('/home/saiant01/Desktop/cat_prodigal_culture_cd-hit.fasta.clstr'))


        #P4J0 vs P4J7
        #graph.load_graph(('/home/saiant01/Desktop/cat_prodigal_cd-hit_p0p7.fasta.clstr'))

        #Local
        #graph.load_graph(('/home/saiant01/Desktop/cat_Sample_P4Jx-Assembly_cd-hit.fa.clstr'))



        #Server
        #graph.load_graph(('/home/saiant01/cat_Sample_PxJy-Assembly_cd-hit.fasta.clstr'))
        #print('Nombre de noeux:'+str(len(graph.nodes)))

        #DataTest
        #graph.load_graph(('/home/saiant01/PycharmProjects/Git/Data/cat_prodigal-cd-hit.fasta.clstr'))
    else:
        print('No graph was given! Use -i to input a file.')



    if find!='':
        #Trouver les clusters pour lequel le gene appartient
        print(graph.find_cluster('Sample_P4J7-FOX-ANA-Assembly.fa_contig-3000007_12'))
    else:
        pass


    #Liste des chemins partant du cluster en faisant au maximum n pas.
    if lenght_lop!=-1 and list_of_paths!='':
        list_of_paths=graph.find_path('Cluster 19164', 10)

        #Coloring
        sequence_path = graph.sequences_in_find_path(list_of_paths)
        #graph.show_path_by_samples(sequence_path)
    else:
        pass

    # Sous-graph
    if sous_graph!='':
        sous_graph = graph.sous_graph(list_of_paths)
    else:
        pass

    #Visualisation
    if cytoscape==True and sous_graph==True:
        sous_graph.cytoscape()
    elif cytoscape==True:
        graph.cytoscape()
    if javacript==True and sous_graph==True:
        sous_graph.graph_javascript()
    elif cytoscape==True:
        graph.graph_javascript()


    #Compare
    if convert_greps_fred!='':
        graph.beta_lactam_file(convert_greps_fred)
    else:
        pass

    if xml!=False and lenght_xml!=-1
        (graph.compare_sequences_excel(xml,lenght_xml))
    else:
        pass


    #Stats
    if stats==True:
        graph.stats_graphs()
    else:
        pass


    #Time
    if time==True:
        print('Time:',graph.function_time(time), '/  H:M:S')