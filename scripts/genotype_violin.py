#!/bin/env python

import pandas as pd
import seaborn as sns
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()


parser.add_argument("--max_cp", type=int, required=False, help="max_cp to plot and filter along", default=-1)
parser.add_argument("--input_table", type=str, required=True, help="Input dataframe")
parser.add_argument("--output_name", type=str, required=True, help="Output filename")
parser.add_argument("--region_name", type=str, required=True, help="Region to plot")
parser.add_argument("--file_type", type=str, required=True, help="pdf|png")
parser.add_argument("--title", type=str, required=True, help="Plot title")
parser.add_argument("--violin", action='store_true', help="Violin plot")
parser.add_argument("--super_pop_only", action='store_true', help="Whether or not to just plot the super population")


args = parser.parse_args()


copy_nums = pd.read_csv(args.input_table, sep='\t',header=0).fillna(0)
filt_copy_nums = copy_nums[copy_nums['name'] == args.region_name]

if len(filt_copy_nums.index) == 0:
	plt.figure(figsize=(12,3))
	plt.ylabel('Super Population')
	plt.xlabel('Copy Number')
	plt.title(args.title)
	plt.text(0.5, 0.5, 'CN Undefined', fontsize=18)
	plt.tight_layout()
else:
	filt_copy_nums['code'] = filt_copy_nums['code'].str.replace('_', ' ')

	filt_copy_nums.loc[:,'super_title'] = ''
	for super_code in filt_copy_nums['code'].unique():
		filt_copy_nums.loc[filt_copy_nums[filt_copy_nums['code'] == super_code].index,'super_title'] = '%s (%d)' % (super_code, len(filt_copy_nums[filt_copy_nums['code'] == super_code].index))

	sns.set_style('whitegrid')

	if args.violin and args.super_pop_only:
		plt.figure(figsize=(12,3))
		sns.violinplot(x="copy_num", y="super_title", color='1.0', data=filt_copy_nums, cut=0.1, inner=None)
		sns.stripplot(x="copy_num", y="super_title", data=filt_copy_nums, jitter=True, hue='super_title', palette='Set1', alpha=0.5)
		plt.xlim(xmin=-0.5, xmax=max(filt_copy_nums['copy_num'])+0.5)
		plt.ylabel('Super Population')
		plt.title(args.title)
		plt.xlabel('Copy Number')
		plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left', ncol=1)
		plt.tight_layout()
		plt.grid(True)
	elif args.super_pop_only and not args.violin:
		sns.stripplot(x="copy_num", y="super_title", data=filt_copy_nums, jitter=True, color='0.0', alpha=0.5)
		plt.xlim(xmin=-0.5, xmax=max(filt_copy_nums['copy_num'])+0.5)
		plt.ylabel('Super Population')
		plt.title(args.title)
		plt.xlabel('Copy Number')
		plt.tight_layout()
		plt.grid(True)
	else:
		fig, ax = plt.subplots(1, 2, figsize=(12,9))
		fig.suptitle(args.title)
		sns.stripplot(x="copy_num", y="super_title", data=filt_copy_nums, jitter=True, color='0.0', alpha=0.5, ax=ax[0])
		ax[0].set_xlim(xmin=-0.5, xmax=max(filt_copy_nums['copy_num'])+0.5)
		ax[0].set_ylabel('Super Population')
		ax[0].set_xlabel('Copy Number')
		sns.stripplot(x="copy_num", y="pop", data=filt_copy_nums, jitter=True, color='0.0', alpha=0.5, ax=ax[1])
		ax[1].set_xlim(xmin=-0.5, xmax=max(filt_copy_nums['copy_num'])+0.5)
		ax[1].set_ylabel('Population')
		ax[1].set_xlabel('Copy Number')
		fig.subplots_adjust(wspace=0.4)
		ax[0].grid(True)
		ax[1].grid(True)

plt.savefig(args.output_name, type=args.file_type)


