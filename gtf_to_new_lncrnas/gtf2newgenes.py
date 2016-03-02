# -*- coding: utf-8 -*-
#!/usr/bin/python
import sys
import os
import argparse
from time import gmtime, strftime
import re



class log():
    def __init__(self, out_file):
        self.out = open(out_file, 'w')
        self.__call__('Log recording start: {} '.format(strftime("%Y-%m-%d %H:%M:%S", gmtime())), 0)

    def __call__(self, msg, level):
        msg_to_out = strftime("%Y-%m-%d %H:%M:%S", gmtime()) + '--> ' + level * 4 * ' ' + msg + '\n'
        sys.stderr.write(msg_to_out)
        self.out.write(msg_to_out)


def main():
    PYTHON_PATH = '/home/lucassilva/virtualenviroment/bin/python'
    GTF_TO_FASTA_PATH = '/usr/bin/gtf_to_fasta'
    LNCSCORE_PATH = '/work/lncrnas/lucassilva/trofoblasto/scripts/lncScore/lncScore.py'



    obj = argparse.ArgumentParser(description='Transform a gtf cufflinks out in a repertory with new lncRNAs\n'
                                              'Obs: This scripts removes transcripts that length less than  200 bp\n'
                                              'Use Example:\n'
                                              "python  path/to/gtf2newgenes.py \n -O path/to/create/the/out_test\n-T path/to/transcripts.gtf\n-G path/to/genome.fa"
                                              )
    obj.add_argument('-G', help='Genome fasta file', required=True)
    obj.add_argument('-T', help=' GTF with Transcripts cufflinks (non-restricted) out', required=True)
    obj.add_argument('-P', help='Number of threads', required=True)
    obj.add_argument('-O', help='Out dir: This dir will be create', required=True)
    parser = obj.parse_args()
    create_out_dir(parser.O)
    log_print = log(parser.O + '/run.log')
    log_print('retirar genes com fpkm 0 ', 5)
    log_print('Converting {} to .FASTA'.format(parser.T), 1)
    gtf_to_fasta(gtf2fasta=GTF_TO_FASTA_PATH,
                 gtf=parser.T,
                 genome=parser.G,
                 out_file=parser.O + '/temp_fasta.fasta',
                 log_print=log_print)
    log_print('Successful conversion'.format(parser.T), 1)

    log_print('Searching for lncRNAs', 1)

    lncscore(python_path=PYTHON_PATH,
             lncscore_path=LNCSCORE_PATH,
             fasta_trans_gtf=parser.O + '/temp_fasta.fasta',
             gtf_ref=parser.T,
             out='lncrnascore_out',
             threads=parser.P,
             hexamer='/work/lncrnas/lucassilva/trofoblasto/scripts/lncScore/dat/Human_Hexamer.tsv',
             train_set='/work/lncrnas/lucassilva/trofoblasto/scripts/lncScore/dat/Human_training.dat',
             log=log_print)

    create_a_gtf_modify(gtf_name=parser.T, saida=parser.O + '/gtf_file_with_predictions.gtf', score_file='lncrnascore_out')
    log_print('Creating a new gtf file with prediction references at {}'.format(parser.O + '/gtf_file_with_predictions.gtf'), 1 )
    os.system('rm {}'.format(parser.O + '/temp_fasta.fasta'))
    log_print('Removing file: {}'.format(parser.O + '/temp_fasta.fasta'), 1 )
    log_print('Analysis successful', 1 )

def create_a_gtf_modify(gtf_name, saida, score_file):
    # print open(score_file).read().split('\n')[1:]
    # print saida
    # print gtf_name

    classification_score = {out_data.split('\t', 1)[0]: out_data.split('\t', 1)[1] for out_data in
                            open(score_file).read().split('\n')[1:] if out_data}
    modify_gtf = []
    for transcript_line in [line for line in open(gtf_name).read().split('\n')[1:]if line]:
        if transcript_line.startswith("#"):
            modify_gtf.append(transcript_line)
        else:
            #print "<"+transcript_line
            chave = re.search('transcript_id \"(\S*)\"', transcript_line).group(1)
            if chave in classification_score:
                transcript_line += ' lncScoreClassification "{}";'.format(classification_score[chave].replace('\t', ':'))
                modify_gtf.append(transcript_line)
            else:
                transcript_line += ' lncScoreClassification "{}";'.format('Not classified')
                modify_gtf.append(transcript_line)

    save_gtf = open(saida, 'w')
    save_gtf.write('\n'.join(modify_gtf))
    save_gtf.close()


def create_out_dir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def gtf_to_fasta(gtf2fasta, gtf, genome, out_file, log_print):
    cmd = '{} {} {} {}'.format(gtf2fasta,
                               gtf,
                               genome,
                               out_file)
    log_print(cmd, 2)

    if os.system(cmd) != 0:
        log_print('Error in fasta conversion...', 0)
        sys.exit(1)
    else:
        temp_to_modify = open(out_file, 'r').read()  #Cleaning transcript name

        temp_to_save = re.sub(r'>\S+\s(\S+).+\n', r">\1\n", temp_to_modify)
        temp_to_save = re.sub(r'[ACTGactg]\n[ACTGactg]', r'', temp_to_save) + '\n'

        temp_to_save = '>'.join([f_line for f_line in temp_to_save.split('>') if len(f_line) > 105])
        temp_to_save = re.sub('\n\n', '\n', temp_to_save)

        # temp_to_save = '>'.join([f_line  if len(f_line)>= 105])

        temp_to_modify = open(out_file, 'w')
        temp_to_modify.write(">" + temp_to_save)

        #print temp_to_save.split('>')[2]


def lncscore(python_path, lncscore_path, fasta_trans_gtf, gtf_ref, out, threads, hexamer, train_set, log):
    cmd = "{} {} -f {} -g {} -o {} -p {} -x {} -t {}".format(python_path,
                                                             lncscore_path,
                                                             fasta_trans_gtf,
                                                             gtf_ref,
                                                             out,
                                                             threads,
                                                             hexamer,
                                                             train_set,
                                                             log)

    log(cmd, 2)
    if os.system(cmd) != 0:
        log('Error in prediction...', 0)
        sys.exit(1)
    else:
        pass


if __name__ == '__main__':
    sys.exit(main())
