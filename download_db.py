#!/usr/bin/env python
# coding: utf-8

# ### PART 1

# In[1]:


####VARIABLES CONFIG 
import statistics
import math
import mysql.connector
from config import *


def get_canonical_transcripts():
    """ 
        Load canonical transcripts from database
    """

    query = "SELECT distinct(FEATURE) from CHIP_DB.TARGET_PANEL"
    mydb = mysql.connector.connect(host=host, charset=charset, user=user, port=port, password=passwd)
    mycursor = mydb.cursor()
    mycursor.execute(query)
    return [record[0] for record in mycursor.fetchall()]

def get_artifacts():
    
    query = "SELECT VAR_ID, GROUP_CONCAT(distinct(PROJECT)) AS PROJECT FROM CHIP_DB.INTERNALLY_IDENTIFIED_NO_DRIVER GROUP BY VAR_ID"
    mydb = mysql.connector.connect(host=host, charset= charset, user=user, port=port, password=passwd)
    mycursor = mydb.cursor()
    mycursor.execute(query)
    return {record[0]: record[1] for record in mycursor.fetchall()}

def get_internally_identified():
    
    query = "SELECT VAR_ID, GROUP_CONCAT(distinct(PROJECT)) AS PROJECT FROM CHIP_DB.INTERNALLY_IDENTIFIED GROUP BY VAR_ID"
    mydb = mysql.connector.connect(host=host, charset= charset, user=user, port=port, password=passwd)
    mycursor = mydb.cursor()
    mycursor.execute(query)
    return {record[0]: record[1] for record in mycursor.fetchall()}

def get_whitelist_ncl():
    
    query = "SELECT CONCAT_WS('_', CHR, POS, REF, VAR) FROM CHIP_DB.WHITELIST_NUCLEOTIDE"
    mydb = mysql.connector.connect(host=host, charset= charset, user=user, port=port, password=passwd)
    mycursor = mydb.cursor()
    mycursor.execute(query)
    return [record[0].decode('ascii') for record in mycursor.fetchall()]

def get_whitelist_aa():

    query = "SELECT FEATURE_PREFIX, AA FROM CHIP_DB.WHITELIST_AMINOACID"
    mydb = mysql.connector.connect(host=host, charset= charset, user=user, port=port, password=passwd)
    mycursor = mydb.cursor()
    mycursor.execute(query)
    return {record[0]: int(record[1]) for record in mycursor.fetchall()}

def get_previously_identified():
    query = "SELECT CONCAT_WS('_', CHR, POS, REF, VAR) FROM CHIP_DB.PREVIOUSLY_IDENTIFIED"
    mydb = mysql.connector.connect(host=host, charset= charset, user=user, port=port, password=passwd)
    mycursor = mydb.cursor()
    mycursor.execute(query)

    return [record[0].decode('ascii') for record in mycursor.fetchall()]



with open("artifacts.txt", 'w') as f:

	f.write(str(get_artifacts()))

with open("internally_identified.txt", 'w') as f:

	f.write(str(get_internally_identified()))

with open("whitelist_ncl.txt", 'w') as f:

	f.write(str(get_whitelist_ncl()))

with open("whitelist_aa.txt", 'w') as f:

	f.write(str(get_whitelist_aa()))

with open("previously_identified.txt", 'w') as f:

	f.write(str(get_previously_identified()))

with open("canonical_transcripts.txt", 'w') as f:

	f.write(str(get_canonical_transcripts()))

