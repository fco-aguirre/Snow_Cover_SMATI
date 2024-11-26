# #!/usr/bin/env python3
# This script orders the anual downloaded files, from the Earthdata, for the correct execution of the main scripts application
# If you like to used the download SMATI location you can left the source and deatination argument as '.'

import argparse
import sys
import os
import shutil

DESC = "This script orders the anual downloaded files, from the Earthdata, for the correct execution of the main scripts application"


def _order(source, dest_1, year, check_name):
    folder_MOD09GA = 'MOD09GA'
    folder_MOD09GQ = 'MOD09GQ'
    folder_MOD35 = 'MOD35'
    #Folder_1 = 'Brunswick_Final'
    #Folder_2 = year
    #Folder_3 = 'Brunswick'
    #archivolec = 'checksums_502266192' # Comment after the test!

    File_ord = 'Order_files'

    if source == '.':
        source_1 = os.getcwd()
        source = source_1
        #source_2 = source_1
        source = os.chdir("..")
        source = os.chdir('./4_Example_Data/Brunswick/' + year)
        source = os.getcwd()
        
        destination = os.chdir("..")
        destination = os.chdir("..")
        destination = os.chdir("..")
        destination = os.chdir("..")
        #path_d = './Outputs/Order_files/' + destination + '/Reflectance_bands'
        destination = os.chdir('./DATA/' + dest_1 + '/Reflectance_bands') # destination by default
        destination = os.getcwd()
        #print(dest_1)

        
    else:
        #destination = os.chdir("..")
        #destination = os.chdir("..")
        #destination = os.chdir("..")
        #destination = os.chdir("..")
        destination = os.chdir(dest_1)
        destination = os.getcwd()
        
        
        
    path_1 = source + '/' + check_name

    path_2 = destination
    
    '''
        #source_1 = os.getcwd()
        source = source
        destination = os.chdir("..")
        destination = os.chdir("..")
        destination = os.chdir('./Outputs/Order_files/' + destination + '/Reflectance_bands') # destination by default
        destination = os.getcwd()
        '''


    
    

    #path_3 = os.path.join(os.path.expanduser('~'), source)

    

    file = open(path_1, 'r')
    linea_file =[]

    narchivos_MOD09GA = []
    narchivos_MOD09GQ = []
    narchivos_MOD35 = []

    cont_file = 0
    #cont_MOD09GQ = 0
    #cont_MOD35 = 0

    while True:
        linea = file.readline()  # lee línea
        a = linea.find('MOD09GA')
        d = linea.find('MOD09GQ')
        e = linea.find('MOD35')
        #print(a)
        #print(linea)

        if not linea:
            break  # Si no hay más se rompe bucle
            #print(linea_h)

        linea_file.append(linea)

        if a >= 0:
            b = len(linea_file[cont_file])

            c = linea_file[cont_file]

            narchivos_MOD09GA.append(c[a:b - 1])

            # print(narchivos[cont])
        elif d >= 0:
            b = len(linea_file[cont_file])

            c = linea_file[cont_file]

            narchivos_MOD09GQ.append(c[d:b - 1])
    
        elif e >= 0:
            b = len(linea_file[cont_file])

            c = linea_file[cont_file]

            narchivos_MOD35.append(c[e:b - 1])

        cont_file = cont_file + 1

    # create sub folders
    year_dir = os.mkdir(path_2 + '/' + year)
    MOD09GA_dir = os.mkdir(path_2 + '/' + year + '/MOD09GA')
    MOD09GQ_dir = os.mkdir(path_2 + '/'+ year+'/MOD09GQ')
    MOD35_dir = os.mkdir(path_2 + '/' + year + '/MOD35')

    #MOD09GA Files
    for k in narchivos_MOD09GA:
        #narchivos_year.append(k)
        #l = l + 1
        #print(k)
        shutil.copy2(source + '/' + k,
                     path_2 + '/' + year + '/MOD09GA')

    #MOD09GQ Files
    for j in narchivos_MOD09GQ:
        #narchivos_year.append(k)
        #l = l + 1
        #print(k)
        shutil.copy2(source + '/' + j, 
                     path_2 + '/' + year + '/MOD09GQ')

    #MOD35 Files
    for l in narchivos_MOD35:
        #narchivos_year.append(k)
        #l = l + 1
        #print(k)
        shutil.copy2(source + '/' + l,
                     path_2 + '/' + year + '/MOD35')


    # write the yearly master files
    with open(path_2 + '/' + year + '/MOD09GA/'+ year +'_MOD09GA', 'w') as f:
        f.writelines("%s\n" % lin for lin in narchivos_MOD09GA)

    f.close()

    with open(path_2 + '/'+ year +'/MOD09GQ/'+ year +'_MOD09GQ', 'w') as f:
        f.writelines("%s\n" % lin for lin in narchivos_MOD09GQ)

    f.close()

    with open(path_2 + '/' + year +'/MOD35/' + year + '_MOD35', 'w') as f:
        f.writelines("%s\n" % lin for lin in narchivos_MOD35)

    f.close()


    #return print(path_1, path_2)


def _main(argv):
    #parser = argparse.ArgumentParser(description = 'Order Files Script')
    parser = argparse.ArgumentParser(prog=argv[0], description=DESC)
    parser.add_argument('-s', '--source', metavar='DIR', help = 'the input folder path', required=True)
    parser.add_argument('-d', '--destination', metavar='DIR', help = 'the destination folder name', required=True)
    parser.add_argument('-y', '--year', help = 'year to order', required=True)
    parser.add_argument('-n', '--check_name', help = 'checksums name', required=True)
    #self.inp_path = parser.add_argument('-inp_p', help = "put your input path")

    args = parser.parse_args(argv[1:])
    '''
    if not os.path.exists(args.destination):
        os.makedirs(args.destination)
    return args.source, args.destination, args.year, args.check_name
    '''
    return _order(args.source, args.destination, args.year, args.check_name)


if __name__ == '__main__':
    try:
        #print('hola')
        #print(sys.argv[0])
        sys.exit(_main(sys.argv))
        #print(sys.exit(_main(sys.argv)))
        #_main()

        #_order(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

        #print(sys.argv)
    except KeyboardInterrupt:
        sys.exit(-1) # tells the program to quit
        #print('hola2')
    
    #print('hola')
