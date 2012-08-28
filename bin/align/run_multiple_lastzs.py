#!/usr/bin/env python
# encoding: utf-8

"""
run_multiple_lastzs.py

Created by Brant C. Faircloth on 07-30-2010.
Copyright (c) 2010 Brant Faircloth and collaborators. 
All rights reserved.

USAGE:  python run_multiple_lastzs.py \
    --configuration=../Simulation/db.conf \
    --output=../Simulation/Data/New/ \
    --probefile=../Simulation/SureSelectProbes.2bit \
    --chromolist='oviAri1' \
    --readlist='ailMel1','myoLuc1'

"""

import pdb
import os
import sys
import oursql
import optparse
import ConfigParser

def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--configuration', dest = 'conf', action='store', \
type='string', default = None, help='The path to the configuration file.', \
metavar='FILE')
    p.add_option('--output', dest = 'output', action='store', \
type='string', default = None, help='The path to the output directory.')
    p.add_option('--probefile', dest = 'probefile', action='store', \
type='string', default = None, help='The path to the probe file.')
    p.add_option('--chromolist', dest = 'chromolist', action='store', \
type='string', default = None, help='The list of organisms with chromosome type genomic sequences.')
    p.add_option('--readlist', dest = 'readlist', action='store', \
type='string', default = None, help='The list of organisms with read type genomic sequences (e.g. uses --huge option).')
    p.add_option('--db', dest = 'db', action='store_true', \
default = False, help='Insert results to exiting database.')
    p.add_option('--fish', dest = 'fish', action='store_true', \
default = False, help='If using fish directories.')
    (options,arg) = p.parse_args()
    if not options.probefile:
        print "You must provide a valid path to the probe file."
        p.print_help()
        sys.exit(2)
    if not os.path.isdir(options.output):
        print "You must provide a valid path to the output directory."
        p.print_help()
        sys.exit(2)
    if options.chromolist:
        options.chromolist = options.chromolist.strip('\n').split(',')
    if options.readlist:
        options.readlist = options.readlist.strip('\n').split(',')
    if options.fish:
        options.fish = "fish/"
    else:
        options.fish = ''
    return options, arg 

def rest_of_steps(cur, g, output_file):
    print "Running sed against {0}".format(output_file)
    new_output_file = output_file + '.clean'
    sed_str = "sed 's/%//g' {0} > {1}".format(output_file, new_output_file)
    os.system(sed_str)
    print "Creating {0} table".format(g)
    mysql_str = '''
        create table {0} (
        id mediumint unsigned not null AUTO_INCREMENT,
        score int unsigned not null,
        name1 varchar(100) not null,
        strand1 varchar(1) not null,
        zstart1 int unsigned not null,
        end1 int unsigned not null,
        length1 smallint unsigned not null,
        name2 text not null,
        strand2 varchar(1) not null,
        zstart2 smallint unsigned not null,
        end2 smallint unsigned not null,
        length2 smallint unsigned not null,
        diff text not null,
        cigar varchar(50) not null,
        identity varchar(8),
        percent_identity float,
        continuity varchar(8),
        percent_continuity float,
        coverage varchar(8),
        percent_coverage float,
        PRIMARY KEY (id),
        INDEX name2 (name2(50))) ENGINE=InnoDB;'''.format(g)
    try:
        cur.execute(mysql_str)
    except oursql.ProgrammingError, e:
        if e == 1050:
            drop_table = "DROP TABLE {0}".format(g)
            cur.execute(drop_table)
            cur.execute(mysql_str)
    print "Inserting data to {0} table".format(g)
    mysql_str2 = '''load data local infile '{0}' 
    into table {1} (score, name1, strand1, zstart1, end1, length1, 
    name2, strand2, zstart2, end2, length2, diff, cigar, identity, 
    percent_identity, continuity, percent_continuity, coverage, 
    percent_coverage);'''.format(new_output_file, g)
    cur.execute(mysql_str2)

def main():
    conf = ConfigParser.ConfigParser()
    options, args = interface()
    if options.conf:
        conf.read(options.conf)
    # build our configuration
    if options.db:
        conn = oursql.connect(
        user=conf.get('Database','USER'),
        passwd=conf.get('Database','PASSWORD'),
        db=conf.get('Database','DATABASE')
        )
        cur = conn.cursor()
    #genomes = ['hg19', 'venter', 'chinese', 'korean', 'panTro2']
    #genomes = [
    #'129S1_SvImJ_Mouse_Genome',
    #'129S5_Mouse_Genome',
    #'C57BL_6N_Mouse_Genome',
    #'CAST_Ei_Mouse_Genome',
    #'NOD_Mouse_Genome',
    #'NZO_Mouse_Genome',
    #'PWK_Ph_Mouse_Genome',
    #'Spretus_Ei_Mouse_Genome',
    #'WSB_Ei_Mouse_Genome',
    #'mm9',
    #'rn4'
    #]
    if options.db:
        try:
            cur.execute('''CREATE TABLE species (
            name varchar(7) NOT NULL,
            description varchar(100) NULL,
            version text NULL,
            PRIMARY KEY (name)) ENGINE=InnoDB''')
        except:
            # need better handling here for tables that exist
            pass
    if options.readlist:
        for g in options.readlist:
            #pdb.set_trace()
            output_file = os.path.abspath(os.path.join(options.output,       
                        "all_probes_v_{0}.lastz".format(g)))
            exc_str = '''/Users/bcf/git/brant/seqcap/Alignment/run_lastz.py \
                        --target=/nfs/data1/genomes/Genomes/{3}{0}/{0}.2bit \
                        --query={1}\
                        --nprocs=6 \
                        --output={2} --huge'''.format(g, options.probefile, output_file, options.fish)
            os.system(exc_str)
            if options.db:
                rest_of_steps(cur, g, output_file)
                try:
                    cur.execute('INSERT INTO species (name) VALUES (?)', (g,))
                except oursql.IntegrityError, e:
                    if e == 1062:
                        pass
    if options.chromolist:
        for g in options.chromolist:
            output_file = os.path.abspath(os.path.join(options.output,       
                        "all_probes_v_{0}.lastz".format(g)))
            exc_str = '''/Users/bcf/git/brant/seqcap/Alignment/run_lastz.py \
                        --target=/nfs/data1/genomes/Genomes/{3}{0}/{0}.2bit \
                        --query={1}\
                        --nprocs=6 \
                        --output={2}'''.format(g, options.probefile, output_file, options.fish)
            os.system(exc_str)
            if options.db:
                rest_of_steps(cur, g, output_file)
                try:
                    cur.execute('INSERT INTO species (name) VALUES (?)', (g,))
                except oursql.IntegrityError, e:
                    if e == 1062:
                        pass
    if options.db:
        conn.commit()
        cur.close()
        conn.close()

if __name__ == '__main__':
    main()
