#!/usr/bin/env python
# coding: utf-8


import statistics
import math
import mysql.connector
import pandas as pd
from config import *


def get_canonical_transcripts():
	with open("canonical_transcripts.txt", 'r') as f:
		return eval(f.read())
def get_artifacts():
        with open("artifacts.txt", 'r') as f:
            return eval(f.read())

def get_internally_identified():
    with open("internally_identified.txt", 'r') as f:
        return eval(f.read())

def get_whitelist_ncl():
    with open("whitelist_ncl.txt", 'r') as f:
        return eval(f.read())

def get_whitelist_aa():
    with open("whitelist_aa.txt", 'r') as f:
        return eval(f.read())

def get_previously_identified():
    with open("previously_identified.txt", 'r') as f:
        return eval(f.read())

def get_longitudinal_pairs():
    data = pd.read_csv("longitudinal_pairs.txt",  header=None, sep='\t', comment='#', na_values='Nothing')
    visit_1 = list(data.loc[:, 1])
    visit_2 = list(data.loc[:, 2].astype(str))
    return visit_1, visit_2

def get_cosmic():
        """
            Create COSMIC table from the file
        """
        res = dict()
        try:
            file = 'COMMON_IN_COSMIC.tsv' ##cambiar con base de datos cosmic
            for line in open(file):
                fields = line.split()
                if len(fields) == 4:
                    table_key = f"{fields[0]}_{fields[1]}_{fields[2]}"
                    res[table_key] = int(fields[3])
            return res
        except Exception as e:
            print(f'ERROR loading COSMIC table : {e}')
            exit(1)

def create_translation_table_from_file():
           res = dict()
           translation_filename = 'conversion_table.tsv'
           try:
               with open(translation_filename, 'r') as translation_file:
                   for line in translation_file:
                       orig, dest = line.strip().split('\t')
                       res[orig] = dest
               return res
           except OSError as e:
               print(f"Unable to open {translation_filename}: {e}", file=sys.stderr)
               exit(1)
 
       
translation_table = create_translation_table_from_file()
       
cosmic_table = get_cosmic()

canonical_transcripts =  get_canonical_transcripts()

artifacts = get_artifacts()

internally_identified = get_internally_identified()

whitelist_ncl = get_whitelist_ncl()
    
whitelist_aa = get_whitelist_aa()

previously_identified = get_previously_identified()

canonical_transcripts = get_canonical_transcripts()

visit_1, visit_2 = get_longitudinal_pairs()

first_line = 1


# In[2]:


def determine_data_type(value):
    """
    The function takes a string input and determines its data type to be either a float, int, or string. 
    """
    try:
        int(value)
        return (int)
    except:
        try:
            float(value)
            return (float)
        except:
            str(value)
            return (str)
    raise NotImplementedError()
    


# In[3]:


#

def get_annotation(vcf_chrom, vcf_pos, vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, csq_header, variante):
        """
            Devuelve la 'anotación principal'. Se considera 'anotacion principal' la correspondiente al transcrito
            canonico y ,en su defecto, a la más deleterea de las anotaciones (una de ellas si hay varias)
            Se incluyen campos adicionales con los id de todas las anotaciones deletereas y sinonimas
        """

        def get_var_data():
            """  Devuelve un dict con los datos del campo INFO. El campo CSQ se modifica para que
                la `key` se corresponda con los alelos del campo ALT
            """

            def unzip_info():
                """
                    Devuelve un dict para todos los items del campo INFO
                """

                def unzip_single_info_item(item):
                    ''' Devuelve un dict para un item del campo INFO '''

                    tmp = item.split('=', 1)
                    if len(tmp) == 1:
                        res = {tmp[0]: True}
                    else:
                        # Si es el campo CSQ, lo desglosamos
                        if tmp[0] == 'CSQ':
                            # Obtenemos las anotaciones de los diferentes transcritos por cada alelo
                            """
                                csq = dict()
                                for tmp_item_annotation in tmp[1].split(','):
                                    item_annotation = dict(zip(csq_header, tmp_item_annotation.split('|')))
                                    allele = item_annotation.get('Allele')
                                    if item_annotation.get('Feature_type') == 'Transcript':
                                        csq.setdefault(allele, []).append(item_annotation)
                                    else:
                                        csq.setdefault(allele, [])
                                res = {tmp[0]: csq}
                            """
                            """
                                csq = dict()
                                for tmp_item_annotation in tmp[1].split(','):
                                    item_annotation = dict(zip(csq_header, tmp_item_annotation.split('|')))
                                    allele = item_annotation.get('Allele ')
                                    if item_annotation.get(' Feature_type ') == 'transcript':
                                        csq.setdefault(allele, []).append(item_annotation)
                                    else:
                                        csq.setdefault(allele, []).append(item_annotation)
                                res = {tmp[0]: csq}
                            """
                            csq = dict()
                            ann = []
                            j=0
                           
                            for tmp_item_annotation in tmp[1].split(','):
                                item_annotation = list(tmp_item_annotation.split('|'))
                                ann.append(item_annotation)
                                i=0
                                for label in csq_header:
                                    csq.setdefault(label, []).append(ann[j][i])
                                    i=i+1
                                j=j+1
                                 
                                
                            res = {'CSQ': csq}
                        else:
                            res = {tmp[0]: tmp[1]}

                    return res

                res = dict()
                {res.update(unzip_single_info_item(item)) for item in vcf_info.split(';')}
                return res

            def match_csq_allele_annotation(csq):
                """
                    Devuelve dict donde la `key` se corresponde con el el campo `alt` del VCF.
                    Modifica las `key` si es necesario, y las deja igual si no lo es

                    NOTA de la especificacion VCF:

                        For simple insertions and deletions in which either the REF or one of the
                        ALT alleles would otherwise be null/empty, the REF and ALT Strings must
                        include the base before the event, unless the event
                        occurs at position 1 on the contig in which case it must include the base after
                        the event; this padding base is not required (although it is permitted)
                        for e.g. complex substitutions or other events where all alleles have at least
                        one base represented in their Strings.
                """

                def is_padded():
                    """
                        Devuelve un booleano indicando si TODOS los alt (el el ref) son padded

                        NOTA: Parece (pero no estoy seguro) que VEP solo anota el alelo trimado (eliminado el padding)
                            si todos los alt tienen padding.
                    """
                    return all(alt_allele[idx_padding] == vcf_ref[idx_padding] for alt_allele in vcf_alt.split(','))

                def update_allele(allele):
                    """
                        Modifica `allele` para que coincida con el campo ALT

                        NOTE: No usamos vcf_alt para generar el alelo a devolver para asegurarnos que
                            lo hacemos correctamnte, de modo que si no hay un match entre el calculado
                            (el que devolvemos) y vcf_alt, salte un error
                    """

                    # En este caso, como todavia no se contempla la causistica del padding en la primero posicion
                    # del cromosoma, el caracter de padding es la primera base de vcf_ref y vcf_alt
                    padding_char = vcf_ref[0]
                    if vcf_ref in ['<', '>'] or vcf_alt in ['<', '>']:
                        # Symbolic allele
                        # Desconozco como aparecería en el CSQ; me imagino que igua
                        # Como este caso no lo controlo. Devuelvo ERROR
                        print("ERROR: Symbolic Allele")
                        exit(1)
                    elif allele == '-':
                        # Si el alelo anotado en CSQ es '-' significa que es una deleccion, y por tanto, se debe
                        # establecer la clave de CSQ con el caracter de padding.
                        new_alt = padding_char
                    else:
                        # Hay que modificar la `key` de las anotaciones CSQ porque tb hay padding
                        if allele == None:
                            new_alt = padding_char
                        else:
                            new_alt = padding_char + str(allele)

                    if vcf_alt != new_alt:
                        # Esta comprobacion tiene sentido porque el VCF esta normalizado y, por tanto solo hay un
                        # elemento en vcf_alt y tiene que coincidir con el anotado
                        print(f'ERROR: CSQ annotation do not match VCF ALT field : '
                              f'probably VCF normalized after annotation.')
                        exit(1)
                    else:
                        return new_alt

                # Hago unas comprobaciones de seguridad (DESCOMENTAR 2 primeros if)
                if (len(vcf_ref) > len(vcf_alt)) and (vcf_ref[0] != vcf_alt[0]):
                    print("Deleccion sin padding")
                    exit(1)
                elif (len(vcf_ref) < len(vcf_alt)) and (vcf_ref[0] != vcf_alt[0]):
                    print("Insercion sin padding")
                    exit(1)
                # elif len(vcf_ref) > 1 and len(vcf_alt) > 1:
                #    print("Complex substution")

                res = dict()
                idx_padding = -1 if vcf_pos == 1 else 0
                # Mientras no lo manejo bien (hay que repasar update_allele)
                if idx_padding == -1:
                    print('ERROR: Padding en posicion 1 del chromosoma')
                    exit(1)
                if is_padded():      
                    
                    for i in range(len(csq['Allele'])):
                        csq['Allele'][i] = update_allele(csq['Allele'][i])
                    res = csq
                else:
                    res = csq

                return res

            annotations = unzip_info()
            annotations['CSQ'] = match_csq_allele_annotation(annotations.get('CSQ'))

            return annotations
        
        def get_canonical(a):
            """
            Devuelve la anotacion canonica de entre los genes a considerar. Si hay varias devuelve la deleterea. Si hay varias deletereas devuelve
            la mas deleterea o, en su defecto, una cualquiera.
            """

            def is_canonical(a, i):
                feature_prefix = a.get('Feature')[i].split('.')[0]
                return feature_prefix in canonical_transcripts

            canonical = dict()
            for i in range(len(a['Feature'])):
                if is_canonical(a, i) and a.get('SYMBOL')[i] not in genes_to_exclude() and is_more_deleterious(a, i, canonical):
                    for key in a:
                        
                        canonical[key] = a[key][i]
                        
                        
            return get_dbns_canonical(canonical)

        def get_dbns_canonical(a):
            i = 0
            canonical = dict()
            if 'VEP_canonical' in a:
                for item in a['VEP_canonical'].split('&'):
                    if item == 'YES':
                        for key in a:
                            try:
                                if a[key].split('&')[i]:
                                    canonical[key] = a[key].split('&')[i]
                            except:
                                canonical[key] = a[key]
                        return canonical
                    i=i+1
            i = 0
            if 'Ensembl_transcriptid' in a:
                    for item in a['Ensembl_transcriptid'].split('&'):
                        if item == a['MANE_SELECT'].split('.')[0]:
                            for key in a:
                                try:
                                    if a[key].split('&')[i]:
                                        canonical[key] = a[key].split('&')[i]
                                except:
                                    canonical[key] = a[key].split('&')[0]
                            return canonical
                        i=i+1
            return a
		
        def get_more_deleterious(a):
            """
                Devuelve la anotacion deleterea con mas impacto de entre los genes a considerr.
                Si hay varias de ellas devuelve una.
            """

            anotacion = dict()
            for i in range(len(a['Feature'])):
                if a.get('SYMBOL')[i] not in genes_to_exclude() and is_more_deleterious(a, i, anotacion):
                    for key in a:
                        anotacion[key] = a[key][i]
            return get_dbns_canonical(anotacion)

        def is_more_deleterious(x_annot, i, y_annot):
            """
                    return TRUE if x is more deletereous than y. Any value is more deleterious than None
            """
            consequences = {
                'HIGH': 40,
                'MODERATE': 30,
                'LOW': 20,
                'MODIFIER': 10
            }

            x = 0 if x_annot is None else consequences.get(x_annot.get('IMPACT', None)[i], 0)
            y = 0 if y_annot is None else consequences.get(y_annot.get('IMPACT', None), 0)

            return x > y
            
            
        #canonical_transcripts = get_canonical_transcripts()
        var_info = get_var_data()
        csq = var_info['CSQ']
        
         # Obtiene la anotacion del transcrito canonico

        canonical = get_canonical(csq)
        if canonical:
            anotacion = canonical
        else:
            # Obtiene la anotacion mas deleterea
            anotacion = get_more_deleterious(csq)

        var_info['CSQ'] = anotacion
        
        return var_info
   

    
def get_mutation(alt_field):
    """
        Returns the mutation.

        Exists the program if >1 alt allele because the VCF must be normalized

    """

    item = alt_field.split(',')

    if len(item) > 1:
        print("ERROR: VCF not normalized. More than 1 alt alleles in the VCF")
        exit(1)
    else:
        return item[0]
        
        
def extract_CSQ_header(line):
    """
    Extraxt CSQ header

    :param line: line from the VCF file
    :return: named-array
    """
    idx_start = line.find('Allele')
    return line[idx_start:-2].split('|')



def genes_to_exclude():
    """
        Listado de genes que se anotan (por corrdenadas y que excluimos). Lo hacemos asi porque no tenemos
        el listado de los genes del panel. Lo idel sería emitir solo los del panel y no exluyendo esto. Pero
        por el momento, para salir del paso, lo hacemos asi
    """

    # TODO: Incluir solo los genes del panel. Listado en .ini (si se aplica o no)
    return ('KMT2B', 'C7orf55-LUC7L2', 'LUC7L2', 'HSFX1', 'LOC107985678', 'SMPX')


# In[4]:


def determine_data_type_of_list(values):
    """
    Function takes a list of strings and determines their data type. 

    """
    data_type_dict={'<class \'int\'>':0, '<class \'float\'>':0, '<class \'str\'>':0}
    for i in range(len(values)):
        data_type = (determine_data_type(values[i]))
        data_type_dict[str(data_type)]=data_type_dict[str(data_type)]+1
    if(data_type_dict['<class \'float\'>']==0 and data_type_dict['<class \'str\'>']==0):
        return int
    elif(data_type_dict['<class \'str\'>']==0 and data_type_dict['<class \'float\'>'] > 0):
        return float
    else:
        return str
    raise NotImplementedError()




# In[5]:


def format_sample_fields(format_field, sample_field):
    """
    Formats the sample fields given the description above.
    """
    final_dict ={}
    list_format = format_field.split(':')
    for item in sample_field:
        output_field_dict = {}
        value_sample_field = []
        value_sample_field = sample_field[item].split(':')        
        for i in range(len(list_format)):
            output_field_dict[list_format[i]]=value_sample_field[i]
        final_dict[item]=output_field_dict
    return final_dict
    raise NotImplementedError()



# In[6]:


def create_dict_from_line(header, line):
    """
    Given the header and a single line, transform them into dictionary as described above. Header and line input are provided in this cell. 
    """
    list_line = line.split('\t')
    intermediate_dict={}
    for item in range(len(header)):
        intermediate_dict[header[item]]=list_line[item]
    a=intermediate_dict['FORMAT']
    del intermediate_dict['FORMAT']
    b={}
    sample_list= header[9:]
    for i in range(len(sample_list)):
        b[sample_list[i]]=intermediate_dict[sample_list[i]]
        del intermediate_dict[sample_list[i]]
    sample_dict = format_sample_fields(a,b)
    intermediate_dict['SAMPLE']=sample_dict
    return intermediate_dict
    raise NotImplementedError()




# In[7]:


def read_vcf_file(filename, outputname):
    """
  
    """
    with open(filename) as f:
        total_lines = sum(1 for line in f )
    with open(filename) as f:
        csq_header=''
        global first_line
        with open(outputname, 'w') as g:
            g.write('[')
        with open('annotated_json_for_pandas.json', 'w') as g:
            g.write('[')
        with open('annotated_filtered.json', 'w') as g:
            g.write('[')  
        n_lines = 0
        for line in f:
            if line.startswith('##INFO=<ID=CSQ'):  #CAMBIAR POR CSQ
                # Extract CSQ header from vcf
                csq_header = extract_CSQ_header(line)
            if(line.startswith("##")==False):
                if(line.startswith('#')):
                    line=line.strip('#').strip('\n')
                    header_list=line.split('\t')
                else:
                    # Obtiene la informacion del locus/var
                    vcf_chrom, vcf_pos, vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, *vcf_samples                         = line.lstrip().split()
                    # NIVEL locus/var
                    variante = get_mutation(vcf_alt)
                    var_id = f"""{vcf_chrom}_{vcf_pos}_{vcf_ref}_{variante}"""
                    annotation = get_annotation(vcf_chrom, vcf_pos, vcf_id, vcf_ref, vcf_alt, vcf_qual, vcf_filter, vcf_info, vcf_format, csq_header, variante)
                    stripped_line=line.strip('\n')
                    single_variant = create_dict_from_line(header_list,stripped_line)
                    single_variant['variantInternalId']= var_id
                    single_variant['variante']= variante
                    single_variant['annotation']=annotation
                    single_variant = vcf2json(single_variant)
                    single_variant = vcf2bff(single_variant, 'all_somatics.json')
                    single_variant = second_filter(first_filter(single_variant))
                    json_for_pandas(single_variant)
                    first_line = 0
                    #variant_list.append(single_variant)
            n_lines = n_lines+1
            if 0.20001> n_lines/total_lines > 0.2:
                print('20% complete')
            if 0.40001> n_lines/total_lines > 0.4:
                print('40% complete')
            if 0.80001> n_lines/total_lines > 0.8:
                print('80% complete') 
        with open(outputname, 'a') as g:
            g.write(']')
        with open('annotated_json_for_pandas.json', 'a') as g:
            g.write(']')
        with open('annotated_filtered.json', 'a') as g:
            g.write(']') 
    #return variant_list
    print(total_lines)





# In[8]:


def extract_info_field(data):
    """
    See description in part 6
    """
    return data['INFO']
    raise NotImplementedError()



# In[9]:


def create_dictionary_of_info_field_values(data):
    """
    We now need to figure out that data types for each of the info fields. Below is the function that first takes the info fields and turns them into a dictionary. The function skips any fields that do not have a value or are missing a value. Also replace \x3b with a comma and \x3d with an equal sign.
    """
    keys = []
    values = []
    for x in data[0].split(';'):
        try:
            val = (x.split('=')[1].replace('\\x3b',',').replace('\\x3d','=')).strip('^\.$')
            values.append(val)
            key = x.split('=')[0].replace('\\x3b',',').replace('\\x3d','=')
            keys.append(key)
        except:
            pass
    
    output = dict(zip(keys,([x] if x else [] for x in values)))
    #try:
        #output['ANN'] =  ##CAMBIAR A CSQ
    #except:
        #pass
    
    return output
    raise NotImplementedError()



# In[10]:


def determine_data_type_of_info_fields(data):
    """
    This function's input is the output from `create_dictionary_of_info_field_values` and uses the previously written function `determine_data_type_of_list` to determine the data type of each of the info fields. The output is a dictionary whose keys are the name of the info fields and values are the data type. 
    """
    data_type_dict={}
    for item in data:
        data_value_list = data[item]
        list_data_type = determine_data_type_of_list(data_value_list)
        data_type_dict[item]=list_data_type
    return data_type_dict    
    raise NotImplementedError()




# In[11]:


def format_data(data, info_field_data_type):
    data_list=[]
    info_field_data = extract_info_field(data) # extract all the info fields
    for i in range(len(data)):
        info_field_list_dict={}
        temp_list=[]
        temp_list.append(info_field_data[i])
        info_field_list = create_dictionary_of_info_field_values(temp_list) # create dictionary from info fields
        info_field_data_type = determine_data_type_of_info_fields(info_field_list) # Determine data type of each info field
        for item in info_field_list:
            try:
                info_field_list[item]=(info_field_data_type[item](info_field_list[item][0]))
                info_field_list_dict[item]=info_field_list[item]
            except IndexError as error:
                continue
        data_dict=data
        data_dict['POS']=int(data_dict['POS'])
        data_dict['QUAL']=data_dict['QUAL']
        data_dict['INFO']=info_field_list_dict
    return data_dict
    raise NotImplementedError()

 

# In[12]:


def save_data_as_json(data, output):
    import json
    with open(output,'a') as file:
        if first_line == 1:
            return file.write(str(json.dumps(data, sort_keys=False, indent=4, separators=(',', ': '))))
        else:
            return file.write(','+str(json.dumps(data, sort_keys=False, indent=4, separators=(',', ': '))))
    raise NotImplementedError()


# In[13]:


def load_data_from_json(filename):
    import json
    '''This function whose input is a filename for a json file. The function uses the filename to read the JSON file in which we saved our final parsed data. '''
    with open(filename) as file:
        return json.load(file)
    raise NotImplementedError()


# ### PART 12

# In[14]:


def vcf2json(data): #CHROM, REF, ALT, POS, 
    #data = read_vcf_file(inp) # read vcf file
    info_field_data = extract_info_field(data) # extract all the info fields
    info_field_list = create_dictionary_of_info_field_values(info_field_data) # create dictionary from info fields
    info_field_data_type = determine_data_type_of_info_fields(info_field_list) # Determine data type of each info field
    data = format_data(data, info_field_data_type) # format the data variable -- from data = read_vcf_file(filename)
    #save_data_as_json(data, out) # save the formatted data
    #data_loaded = load_data_from_json(filename) # load saved data
    #final_variant_list=[]
    #for i in range(len(data_loaded)):
    #   if (data_loaded[i]['CHROM']==CHROM and data_loaded[i]['REF']==REF and data_loaded[i]['ALT']==ALT and data_loaded[i]['POS']==POS):
    #        final_variant_list.append(data_loaded[i])
    #    else:
    #        continue
    return data



# In[15]:


def serialize_caseLevelData(data):
    """
   Beacon format caseleveldata field from basic json.
    """
    
    zygosity_dict = [('0/1','GENO_0000458'),
        ('0|1','GENO_0000458'),
        ('1/0','GENO_0000458'),
        ('1|0','GENO_0000458'),
        ('1/1','GENO_0000136'),
        ('1|1','GENO_0000136')]

    def get_sample_name(biosampleId):
       """
       If translation_filename is supplied, get ID from the VCF and return the translated ID, y not supplied, returns
       the original ID
       """
       return translation_table.get(biosampleId, biosampleId) 

    CLD = []
    for individual in data['SAMPLE']:
    
        try:
            if data['SAMPLE'][individual]['GT'] != './.' and data['SAMPLE'][individual]['GT'] != '././.': #and data['SAMPLE'][individual]['AF'] > VAF_THRESHOLD quitado porque no pasa ni un sample
                sample = {}
                
                sample['biosampleId']= get_sample_name(individual)
                sample['zygosity']= {}
                sample['zygosity']['label']= data['SAMPLE'][individual]['GT']
                for key, val in zygosity_dict:
                    if key == data['SAMPLE'][individual]['GT']:
                        sample['zygosity']['id'] = 'GENO:' + val      
                CLD.append(sample)
                for info in data['SAMPLE'][individual]:
                    if info != 'GT':
                        sample[info] = data['SAMPLE'][individual][info]
        except:
            pass
    return CLD
    raise NotImplementedError()
    

# In[16]:


def serialize_frecuencyInPopulation(data): #Solo para una sola notacion de frecuencias 
    frecuencies = dict()
    frecuencies_list = []
    frec_f = []
    temp = dict()
    source = [
        ('gnomAD_exomes','The Genome Aggregation Database (gnomAD exomes)'),
        ('gnomAD_genomes','The Genome Aggregation Database (gnomAD genomes)'),
        ('gnomAD','The Genome Aggregation Database (gnomAD)'),
        ('1000Gp3','The 1000 Genomes Project Phase 3'),
        ('ExAC','The Exome Aggregation Consortium (ExAC)')
    ]
    source_ref = [
        ('gnomAD_exomes','https://gnomad.broadinstitute.org'),
        ('gnomAD_genomes','https://gnomad.broadinstitute.org'),
        ('gnomAD','https://gnomad.broadinstitute.org'),
        ('1000Gp3', 'https://www.internationalgenome.org'),
        ('ExAC', 'https://gnomad.broadinstitute.org')
    ]
    version = [     #RePASAR PORQUE NO SERA LA MISMA VERSION 
        ('gnomAD_exomes', 'Extracted from dbNSFP4.3a'),
        ('gnomAD_genomes', 'Extracted from dbNSFP4.3a'),
        ('gnomAD', 'Extracted from dbNSFP4.3a'),
        ('1000Gp3', 'Extracted from dbNSFP4.3a'),
        ('ExAC', 'Extracted from dbNSFP4.3a')
    ]

    for pop in ['AFR', 'AMR', 'EAS', 'FIN', 'NFE', 'SAS']:
        str_pop = pop + '_AF';
        # For whatever reason freq values are duplicated in some pops (to do: we should check if they're ALWAYS equal)
        if str_pop in data['annotation']['CSQ'] and data['annotation']['CSQ'][str_pop] !='':
            allele_freq = data['annotation']['CSQ'][str_pop]
            tmp = dict()

            tmp['population']= pop
            tmp['alleleFrequency'] = allele_freq
            frecuencies_list.append(tmp)
    if frecuencies_list != []:
        frecuencies['frecuencies']= frecuencies_list
        frec_f.append(frecuencies)
    for key, label in source:
        frecuencies = dict()
        tmp = dict()
        frecuencies_list=[]
        for pop in ['_AFR', '_AMR', '_EAS', '_FIN', '_NFE', '_SAS', '']:
            tmp = dict()
            frecuencies = dict()
            str_pop = key + pop + '_AF';
            # For whatever reason freq values are duplicated in some pops (to do: we should check if they're ALWAYS equal)
            if str_pop in data['annotation']['CSQ'] and data['annotation']['CSQ'][str_pop] !='':
                    allele_freq = data['annotation']['CSQ'][str_pop]
                    tmp = dict()
               
                    tmp['population']= pop[1:] if pop != '' else 'All population'
                    tmp['alleleFrequency'] = allele_freq
                    frecuencies_list.append(tmp)
                    frecuencies['source'] = label
                    frecuencies['sourceReference'] = [label2 for key2, label2 in source_ref if key2==key][0]
                    frecuencies['version'] = [label2 for key2, label2 in version if key2==key][0]
                    frecuencies['frecuencies'] = frecuencies_list
        if frecuencies != {}:
            frec_f.append(frecuencies)
    return frec_f           


def serialize_variation(data):
    
    def guess_variant_type(REF, ALT):
        if len(REF) == len(ALT):
            type = 'SNP'
        else:
            type = 'INDEL'
        return type
    
    variation = dict()
    location = dict()
    interval = dict()
    start = dict()
    end = dict()
    
    
    variation['alternateBases'] = data['ALT']
    variation['referenceBases'] = data['REF']
    if 'VARIANT_TYPE' in data['annotation']:
        variation['variantType'] = data['annotation']['VARIANT_TYPE']
    else:
        variation['variantType'] = guess_variant_type(data['REF'], data['ALT'])
    location['sequence_id'] = "HGSVid:" + str(data['CHROM']) + ':g.' + str(data['POS']) + data['REF'] + '>' + data['ALT']
    location['type'] = 'SequenceLocation'
    interval['type'] = 'SequenceInterval'
    start['type'] = 'Number'
    start['value'] = int(data['POS']) - 1
    end['type'] = 'Number'
    end['value'] = data['POS']
    interval['start'] = start
    interval['end'] = end
    location['interval'] = interval
    variation['location'] = location
    return variation


# In[18]:


def serialize_variantQuality(data):
    quality = dict()
    quality['QUAL'] = data['QUAL']
    quality['FILTER'] = data['FILTER']
    return quality


# In[19]:


def serialize_molecularAttributes(data):
    ma = dict()
    try:
        ma['aminoacidChanges'] = data['annotation']['CSQ']['HGVSp']
    except:
        ma['aminoacidChanges'] = ''
    try:
        ma['geneIds'] = data['annotation']['CSQ']['SYMBOL']
    except:
        ma['geneIds'] = ''
    try:
        ma['genomicFeatures'] = data['annotation']['CSQ']['Feature_type']
    except:
        ma['genomicFeatures'] = ''
    try:
        ma['molecularEffects'] = data['annotation']['CSQ']['IMPACT']
    except:
        ma['molecularEffects'] = ''
    return ma
    


# In[20]:


def serialize_identifiers(data):
    
    #modificar para nuestros datos        
    array_id = [('proteinHGVSIds','Ensembl_proteinid'),
                ('transcriptHGVSIds','_Ensembl_transcriptid')]
    array_alternative_id = [ 
        ('ClinVar','dbNSFP_clinvar_id'),
        ('dbSNP','dbNSFP_rs_dbSNP151')]
    
    identifiers = dict()
    
    if ( 'clinvar_hgvs' in data['annotation']['CSQ'] and data['annotation']['CSQ']['clinvar_hgvs'] != ''):
        identifiers['genomicHGVSId'] = data['annotation']['CSQ']['clinvar_hgvs']
          
    elif ( 'CLINVAR_CLNHGVS' in data['annotation']['CSQ'] and data['annotation']['CSQ']['CLINVAR_CLNHGVS'] != '' ):
        identifiers['genomicHGVSId'] = data['annotation']['CSQ']['CLINVAR_CLNHGVS']
    
    else:
        tmp_id = ':g.' + str(data['POS']) + data['REF'] + '>' + data['ALT']
    
        if 'Ensembl_geneid' in data['annotation']['CSQ'] and data['annotation']['CSQ']['Ensembl_geneid'] != '':
            geneid = data['annotation']['CSQ']['Ensembl_geneid'].split(',')[0]
            geneid = geneid.split('&')[0]
            identifiers['genomicHGVSId'] = geneid + tmp_id
        else:
            identifiers['genomicHGVSId'] = str(data['CHROM']) + tmp_id
            

    # **** clinvarVariantId
    for key, val in array_alternative_id:
        if key=='Clinvar' and val in data['annotation']['CSQ']:
            identifiers['clinvarVariantId'] = key + ':' + data['annotation']['CSQ'][val]
    
    # **** 'proteinHGVSIds' 'transcriptHGVSIds'
    for key, val  in array_id:
        if key =='proteinHGVSIds' or key == 'transcriptHGVSIds':
            if val in data['annotation']['CSQ'] and data['annotation']['CSQ'][val] != '':
                identifiers[key] = data['annotation']['CSQ'][val].split('&')[0]
    
    return identifiers


# In[21]:


import statistics



def serialize_allAnnotations(data, canonical_transcripts, artifacts, internally_identified, whitelist_ncl, whitelist_aa, previously_identified):
    

    def is_previously_identified(var_id):
        return var_id in previously_identified

    
    def is_artifact(var_id):
        return var_id in artifacts
    
    
    def is_whitelist_ncl(var_id):
        return var_id in whitelist_ncl



    def is_whitelist_aa(anotacion):

        def get_aa_position():
            return anotacion.get('Protein_position', '').split('-')[0]
        feature_prefix = anotacion.get('Feature', '').split('.')[0]
        return whitelist_aa.get(feature_prefix) == get_aa_position()

    
    def common_in_cosmic(cosmic_match):
        l = []
        if cosmic_match:
            for i in cosmic_match:
                l.append(i[1])
            return max(l)
    
        
        else:
            return 0
        
    def get_cosmic_match(ids, vcf_ref, variante, cosmic_table):
        
        return [(cosmic_id, cosmic_table.get(f"{cosmic_id}_{vcf_ref}_{variante}"))
                    for cosmic_id in ids if cosmic_table.get(f"{cosmic_id}_{vcf_ref}_{variante}")]
   
        
    def get_existing_variation(anotacion):
        """
            Obtiene los IDs en las BD externas
        """

        def get_ids(id_chunk, mask_chunk):
            # mask = iter(mask_chunk.split('&'))
            mask = mask_chunk.split('&')
            ids = id_chunk.split('&')

            return [ident for ident, flag in zip(ids, mask) if flag == '1']
        
        somatic = get_ids(anotacion['CSQ'].get('Existing_variation'), anotacion['CSQ'].get('SOMATIC')) if 'SOMATIC' in anotacion['CSQ'] else []
        pheno = get_ids(anotacion['CSQ'].get('Existing_variation'), anotacion['CSQ'].get('PHENO')) if 'PHENO' in anotacion['CSQ'] else []

        return somatic, pheno
    
    def predictor_splitting(prediction_chunk):
            """
                Entre la descripcion y el score de la prediccion (p.e SIFT y PolyPhen)
            """
            tmp = prediction_chunk.split('(')

            if prediction_chunk:
                try:
                    desc = tmp[0]
                    score = tmp[1].split(')')[0]
                except IndexError:
                    desc = prediction_chunk
                    score = ''
            else:
                desc = ''
                score = ''

            return desc, score
    def is_vaf_germ(vaf, nocall_as_germ=False):
        """

        Devuelve un boolean indicando si la VAF se considera no compatible con mutacion somatica. El parametro
        `nocall_as_germ` es True si condideramos los nocall como germinal.

        """

        if vaf is None or vaf == '.':
            res = True if nocall_as_germ else False
        else:
            vaf_casted = float(vaf)
            res = (0.45 < vaf_casted < 0.55) or vaf_casted > 0.85

        return res


    def get_samples_with_vaf_like_germinal(samples):  
        return [sample['biosampleId'] for sample in samples if is_vaf_germ(sample['AF'], nocall_as_germ=False)]


    def longitudinal(samples):
        sample_list = []
        visit_1_2 = 0
        only_visit_1 = 0
        only_visit_2 = 0
        for sample in samples:
            sample_list.append(str(sample['biosampleId']))
        for sample in sample_list:
            if sample in visit_1:
                if visit_2[visit_1.index(sample)] in sample_list:
                    visit_1_2 = visit_1_2 + 1 
                else:
                    only_visit_1 = only_visit_1 + 1
            if sample in visit_2:
                if visit_1[visit_2.index(sample)] not in sample_list:
                    only_visit_2 = only_visit_2 + 1
        return [visit_1_2, only_visit_1, only_visit_2]

    
    is_whitelisted = is_whitelist_ncl(data['variantInternalId']) or is_whitelist_aa(data['annotation'])
    samples_with_mutation = serialize_caseLevelData(data)
    vaf_array = [float(i['AF']) for i in samples_with_mutation if 'AF' in i and i['AF'] != '.']
    samples_with_vaf_like_germinal = get_samples_with_vaf_like_germinal(samples_with_mutation)
    mean_vaf = statistics.mean(vaf_array) if vaf_array != [] else 'NA' 
    stdev_vaf = statistics.stdev(vaf_array) if vaf_array != [] and len(vaf_array)>1 else 'NA' 
    max_vaf = max(vaf_array) if vaf_array != [] else 'NA'
    sift_desc, sift_score = predictor_splitting(data['annotation']['CSQ']['SIFT']) if 'SIFT' in data['annotation']['CSQ'] else ['', '']
    polyphen_desc, polyphen_score = predictor_splitting(data['annotation']['CSQ']['PolyPhen']) if 'PolyPhen' in data['annotation']['CSQ'] else ['', '']
    somatic_ids, pheno_ids = get_existing_variation(data['annotation'])
    cosmic_match =  get_cosmic_match(somatic_ids, data['REF'], data['variante'], cosmic_table)
    cosmic_n_match = common_in_cosmic(cosmic_match)
    longi = longitudinal(samples_with_mutation)
    exclude = ['SIFT', 'SOMATIC', 'PHENO', 'PolyPhen', 'AF', 'EAS_AF', 'AMR_AF', '*_AF']
    all_annotations = dict()
    #all_annotations = data['annotation']
    all_annotations['SIFT_DESC'] = sift_desc
    all_annotations['SIFT_SCORE'] =  sift_score
    all_annotations['POLYPHEN_DESC'] = polyphen_desc
    all_annotations['POLYPHEN_SCORE'] = polyphen_score
    all_annotations['MEAN_VAF'] = mean_vaf
    all_annotations['STDEV_VAF'] = stdev_vaf
    all_annotations['MAX_VAF'] = max_vaf
    all_annotations['FRACTION_SAMPLES'] = float(len(samples_with_mutation)/len(data['SAMPLE']))
    all_annotations['NUM_SAMPLES_WITH_MUTATION'] = len(samples_with_mutation)
    all_annotations['NUM_SAMPLES_WITH_VAF_LIKE_GERMINAL'] = len(samples_with_vaf_like_germinal)
    all_annotations['ARTIFACTS'] = is_artifact(data['variantInternalId']) 
    all_annotations['INTERNALLY_IDENTIFIED'] = internally_identified.get(data['variantInternalId'])
    all_annotations['WHITELIST'] = is_whitelisted 
    all_annotations['PREVIOUSLY_IDENTIFIED'] = is_previously_identified(data['variantInternalId'])
    all_annotations['COSMIC_MATCH'] = [].append([i[0] for i in cosmic_match])
    all_annotations['COSMIC_N_MATCH'] = cosmic_n_match
    all_annotations['SOMATIC'] = [].append(somatic_ids) if somatic_ids != '' else None
    all_annotations['PHENO'] = [].append(pheno_ids)  if pheno_ids != '' else None
    all_annotations['LONGITUDINAL_BOTH'] = longi[0]
    all_annotations['LONGITUDINAL_ONLY_VISIT_1'] = longi[1]
    all_annotations['LONGITUDINAL_ONLY_VISIT_2'] = longi[2]
    for key in data['annotation']:
        if key not in exclude:
            if key == 'CSQ':
                for key_CSQ in data['annotation'][key]:
                    if key_CSQ not in exclude:
                        all_annotations[key_CSQ] = data['annotation'][key][key_CSQ]
            else:
                all_annotations[key] = data['annotation'][key]
        
    """


        ,
        'HEURISTIC': 'Y' if heuristic_filter(cosmic_match) else 'N',


        
        'VARSOME': create_varsome_link()
    """
    return all_annotations



# In[23]:


def vcf2bff(json, out):
    variant = dict()
    variant['variantInternalId'] = json['variantInternalId']
    variant['identifiers'] = serialize_identifiers(json)
    variant['variation'] = serialize_variation(json)
    variant['variantQuality'] = serialize_variantQuality(json)
    variant['frecuencyInPopulation'] = serialize_frecuencyInPopulation(json)
    variant['molecularAttributes'] = serialize_molecularAttributes(json)
    variant['caseLevelData'] = serialize_caseLevelData(json)
    variant['_allAnnotations'] = serialize_allAnnotations(json, canonical_transcripts, artifacts, internally_identified, whitelist_ncl, whitelist_aa, previously_identified)
    #variant['_allAnnotations'] = json['annotation']
    #list_variant.append(variant)

    #save_data_as_json(variant, out)
    return(variant)
    


# In[ ]:



def is_in_vaf_threshold_exception_gene(anotacion):
    """
        Listado de genes para los que no se aplica el threshold de la vaf

        #TODO : Listado de genes en el fichero .ini
    """
    return anotacion.get('SYMBOL') in ('DNMT3A', 'TET2', 'ASXL1', 'JAK2', 'TP53', 'PPM1D', 'KDM6A',
                                       'SF3B1', 'SRSF2', 'GNAS', 'GNB1', 'CBL', 'KRAS')


    
def seems_somatic(anotation):
    """
    Una variante se considera somatica si:
        - se detecta en < MAX_NUM_SAMPLES_WITH_VAR muestras
        - Está en menos de MAX_NUM_SAMPLES_GERMINAL_VAF con una VAF aparentemente germinal
        - No está en la población con MAF > MAX_MAF_THRESHOLD

    """

    def is_in_population():
        """ Comprueba si la variante está en la población por encima de una MAF (MAX_MAF_THRESHOLD). """

        # Se produce si se ejecuta VEP con --max_af
        max_af = anotation.get('MAX_AF', -1)
        if max_af == -1:
            print('ERROR: VEP executed witout max_af')
            exit(1)
        elif not max_af:
            return False
        else:
            return float(max_af) > MAX_MAF_THRESHOLD

    ok_samples_with_mutation = anotation['FRACTION_SAMPLES'] < TOTAL_SAMPLES_FRACTION
    ok_samples_with_vaf_like_germinal = anotation['NUM_SAMPLES_WITH_VAF_LIKE_GERMINAL'] < MAX_NUM_SAMPLES_GERMINAL_VAF


    return ok_samples_with_mutation and ok_samples_with_vaf_like_germinal and not is_in_population() 

def first_filter(ej):
    ej['_allAnnotations']['FILTER_1'] = 0
    i=0
    list_filter = ej['variantQuality']['FILTER'].split(';')
    if 'PASS' in list_filter:
        if ej['_allAnnotations']['MAX_VAF'] != 'NA': #si max_vaf es nulo es que no hay samples con vaf mayor al 1%
            ej['_allAnnotations']['FILTER_1'] = 1
        return ej
    try:
        list_filter.remove('strand_bias')

    except:
        i=i+1

    try:
        list_filter.remove('clustered_events')
    except:
        i=i+1
    if list_filter == [] and i != 2:
        if ej['_allAnnotations']['MAX_VAF'] != 'NA': #si max_af es nulo es que no hay samples
            ej['_allAnnotations']['FILTER_1'] = 1
        return ej
        
    return ej

def is_deleterious(a):
 
       def has_impact(annotation):
           '''
               Comprueba si el impacto es HIGH or MODERATE
           '''
 
           impact = annotation['IMPACT']
 
           if not impact:
               print('ERROR: No hay campo IMPACT')
               exit(1)
           else:
               res = impact in ('HIGH', 'MODERATE')
 
           return res
 
       return has_impact(a)


def second_filter(ej):
    ej['_allAnnotations']['FILTER_2'] = 0
    if ej['_allAnnotations']['FILTER_1'] == 1:
        if ej['_allAnnotations']['WHITELIST'] == True and is_deleterious(ej['_allAnnotations']):
            ej['_allAnnotations']['FILTER_2'] = 1

        else:
            if seems_somatic(ej['_allAnnotations']):
                ej['_allAnnotations']['FILTER_2'] = 1
    if ej['_allAnnotations']['FILTER_2'] == 1:
        save_data_as_json(ej, 'annotated_filtered.json')
    return ej

            


# In[ ]:


def json_for_pandas(data):

    variant = dict()
    variant['variantInternalId'] = data['variantInternalId']
    for item in data['_allAnnotations']:
        variant[item] = data['_allAnnotations'][item]


    save_data_as_json(variant, 'annotated_json_for_pandas.json')

read_vcf_file('all_HF_1.vcf', 'annotated_all_somatics.json')

