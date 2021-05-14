#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pandas as pd
from bottle import route, run, template
import requests
    
    
ge = pd.read_csv('gene.expression.csv')
ge['Gene'] = ge.apply(lambda row: f'<a href="/gene/{row.Gene}"">{row.Gene}</a>', axis=1)

te = pd.read_csv('transcript.expression.csv')
te = te.rename(columns={'Gene_ID': 'Gene ID', 'Transcript_ID': 'Transcript ID'})

dc = pd.read_csv('transcript.coordinates.csv')
dc = dc.rename(columns={'Gene_ID': 'Gene ID', 'Transcript_ID': 'Transcript ID'})


def get_content(file):
    with open(file) as f:
        return f.read().strip()
    
    
def get_gene_expressions(n=0):
    de = ge.head(n) if n else ge
    text = de.to_html(classes="table table-striped table-bordered", table_id="expression_table",
                      index=False, escape=False, justify='center')
    return text


def get_table(dx):
    dx = dx.drop(columns=['Gene ID', 'Gene'])
    dx['Transcript ID'] = [f'<a href="/transcript/{t}">{t}</a>' for t in dx['Transcript ID']]
    if dx.shape[0] > 1:
        for column in dx.columns:
            maximum = dx[column].max()
            dx[column] = dx[column].replace(maximum, f'<span class="maximum">{maximum}</span>')
    table = dx.to_html(classes="table table-striped table-bordered", table_id="expression_table",
                       index=False, justify='center', escape=False)
    return table


def get_sequence(chrom, start, end, strand, utr5start, utr5end, utr3start, utr3end):
    server = "http://grch37.rest.ensembl.org"
    # start, end = (utr3start, utr3end) if strand == '+' else (utr3start, utr5end)
    ext = f"/sequence/region/human/{chrom.replace('chr', '')}:{start}..{end}:{1 if strand == '+' else -1}?"
    r = requests.get(server + ext, headers={"Content-Type": "text/plain"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    seq = r.text
    length = len(seq)
    utr5_length, utr3_length = utr5end - utr5start + 1, utr3end - utr3start + 1
    if strand == '+':
        start_codon_start_pos, stop_codon_start_pos = utr5end + 2, utr3start - 2
        start_codon_start, stop_codon_start = utr5_length + 1, length - utr3_length - 2
    else:
        start_codon_start_pos, stop_codon_start_pos = utr5start, utr3end + 4
        start_codon_start, stop_codon_start = utr5_length - 4, length - utr3_length + 3
    sequence = [
        ('utr5', seq[:start_codon_start]),
        ('start_codon', seq[start_codon_start: start_codon_start + 3]),
        ('cds', seq[start_codon_start + 3: stop_codon_start]),
        ('stop_codon', seq[stop_codon_start: stop_codon_start + 3]),
        ('utr3', seq[stop_codon_start + 3:])
    ]
    a = (f'<a class="btn btn-outline-success btn-sm" '
         f'href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={chrom}%3A{start}%2D{end}" '
         f'target="_blank">{chrom}:{start:,}-{end:,}<i class="fas fa-external-link-alt pl-2"></i></a>')
    sequence = '\n'.join([f'<span class="{c}">{t}</span>' for c, t in sequence])
    return sequence, a


def get_transcript_info(gene='', transcript=''):
    dg = te[te['Gene'] == gene].copy() if gene else te[te['Transcript ID'] == transcript].copy()
    if gene:
        transcript = dg.iloc[0]['Transcript ID']
    table = get_table(dg)
    row = dc[dc['Transcript ID'] == transcript].iloc[0]
    gene, gid,  = row['Gene'], row['Gene ID']
    title = f'{gid} ({gene})'
    sequence, a = get_sequence(row['Chrom'], row['Start'], row['End'], row['strand'],
                               row['utr5_start'], row['utr5_end'], row['utr3_start'], row['utr3_end'])
    sequence_header = f'{gid}|{transcript}|{gene} {a}'
    return {'table_title': title, 'table': table,
            'sequence_header': sequence_header, 'sequence': sequence}


@route('/')
@route('/<number>')
def home(number=0):
    table = get_gene_expressions(n=int(number))
    return template('index.html', table=table)


@route('/gene/<gene>')
def gene_expression(gene):
    data = get_transcript_info(gene=gene)
    return template('transcript.html', **data)


@route('/transcript/<transcript>')
def transcript_expression(transcript):
    data = get_transcript_info(transcript=transcript)
    return template('transcript.html', **data)


if __name__ == '__main__':
    run(host='localhost', port=8080, debug=True)
