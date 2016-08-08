import sys

def make_dict_famille(allRG_file, blastp_file,):
    f1=open(allRG_file,'r')
    line=f1.readline()
    line=f1.readline()
    dict={}
    while line!='':
        accession=line.split('\t')[0]
        name=line.split('\t')[6]
        dict[accession]=name
        line=f1.readline()
    f1.close()
    f2=open(blastp_file,'r')
    line=f2.readline()
    while line!='':
        code = (line.split('\t')[1]).split('|')[0]

        #print(line.replace(line.split('\t')[1],dict[code]+'|'+line.split('\t')[1]))
        line=f2.readline()
    f2.close()

if __name__ == '__main__':
    make_dict_famille('/home/saiant01/All-RG_2014-05-23.tsv', '/home/saiant01/Sample_P4Jx-Assembly.blastp')
    #make_dict_famille(sys.argv[1],sys.argv[2])

