#!/bin/env python

import pandas as pd
import seaborn as sns
import argparse
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()


parser.add_argument(
    "--max_cp",
    type=int,
    required=False,
    help="max_cp to plot and filter along",
    default=-1,
)
parser.add_argument("--input_table", type=str, required=True, help="Input dataframe")
parser.add_argument("--output_name", type=str, required=True, help="Output filename")
parser.add_argument("--region_name", type=str, required=True, help="Region to plot")
parser.add_argument("--file_type", type=str, required=True, help="pdf|png")
parser.add_argument("--title", type=str, required=True, help="Plot title")
parser.add_argument("--violin", action="store_true", help="Violin plot")
parser.add_argument(
    "--super_pop_only",
    action="store_true",
    help="Whether or not to just plot the super population",
)


args = parser.parse_args()


copy_nums = pd.read_csv(args.input_table, sep="\t", header=0).fillna(0)
filt_copy_nums = copy_nums[copy_nums["name"] == args.region_name]


if len(filt_copy_nums.index) == 0:
    plt.figure(figsize=(12, 3))
    plt.ylabel("Super Population")
    plt.xlabel("Copy Number")
    plt.title(args.title)
    plt.text(0.5, 0.5, "CN Undefined", fontsize=18)
    plt.tight_layout()
else:
    filt_copy_nums["code"] = filt_copy_nums["code"].str.replace("_", " ")

    color_list = [
        "tab:blue",
        "tab:orange",
        "tab:green",
        "tab:red",
        "tab:purple",
        "tab:pink",
        "tab:gray",
        "tab:cyan",
        "tab:olive",
        "black",
        "gold",
    ]

    superpop_color_dict = {
        "AFR": "gold",
        "EUR": "tab:blue",
        "SAS": "tab:purple",
        "EAS": "tab:green",
        "AMR": "tab:red",
        "OCN": "tab:cyan",
        "EA": "tab:green",
        "MEA": "darkolivegreen",
        "SIB": "firebrick",
        "WEA": "tab:blue",
        "SA": "tab:purple",
        "fa": "darkgreen",
        "mo": "orange",
        "s1": "darkviolet",
        "p1": "teal",
        "ARC": "tab:brown",
        "D&N": "coral",
        "BONOBO": "saddlebrown",
        "CH17": "darkred",
        "CHIMPANZEE": "black",
        "GORILLA": "darkgrey",
        "ORANGUTAN": "darkorange",
        "VMRC59": "darkkhaki",
    }

    color_dict = {}
    extra = 0
    for superpop in filt_copy_nums["code"].unique():
        if superpop in superpop_color_dict:
            superpop_label = "%s (%d)" % (
                superpop,
                len(filt_copy_nums[filt_copy_nums["code"] == superpop].index),
            )
            color_dict[superpop_label] = superpop_color_dict[superpop]
        else:
            superpop_label = "%s (%d)" % (
                superpop,
                len(filt_copy_nums[filt_copy_nums["code"] == superpop].index),
            )
            color_dict[superpop_label] = color_list[extra]
            extra += 1

    filt_copy_nums.loc[:, "super_title"] = ""
    for super_code in filt_copy_nums["code"].unique():
        filt_copy_nums.loc[
            filt_copy_nums[filt_copy_nums["code"] == super_code].index, "super_title"
        ] = "%s (%d)" % (
            super_code,
            len(filt_copy_nums[filt_copy_nums["code"] == super_code].index),
        )

    sns.set_style("whitegrid")

    filt_copy_nums = filt_copy_nums.sort_values(
        by="code",
        key=lambda x: x.map({k: v for v, k in enumerate(superpop_color_dict)}),
    ).reset_index(drop=True)

    if args.violin and args.super_pop_only:
        plt.figure(figsize=(12, 3))
        sns.violinplot(
            x="copy_num",
            y="super_title",
            color="1.0",
            data=filt_copy_nums,
            cut=0.1,
            inner=None,
        )
        sns.stripplot(
            x="copy_num",
            y="super_title",
            data=filt_copy_nums,
            jitter=True,
            hue="super_title",
            palette=color_dict,
            alpha=0.5,
        )
        plt.xticks(ticks=range(0, int(max(filt_copy_nums["copy_num"]))+2))
        plt.xlim(xmin=-0.5, xmax=max(filt_copy_nums["copy_num"]) + 0.5)
        plt.ylabel("Super Population")
        plt.title(args.title)
        plt.xlabel("Copy Number")
        plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", ncol=1)
        plt.tight_layout()
        plt.grid(True)
    elif args.super_pop_only and not args.violin:
        sns.stripplot(
            x="copy_num",
            y="super_title",
            data=filt_copy_nums,
            jitter=True,
            color="0.0",
            alpha=0.5,
        )
        plt.xticks(ticks=range(0, int(max(filt_copy_nums["copy_num"]))+2))
        plt.xlim(xmin=-0.5, xmax=max(filt_copy_nums["copy_num"]) + 0.5)
        plt.ylabel("Super Population")
        plt.title(args.title)
        plt.xlabel("Copy Number")
        plt.tight_layout()
        plt.grid(True)
    else:
        fig, ax = plt.subplots(1, 2, figsize=(12, 9))
        fig.suptitle(args.title)
        sns.stripplot(
            x="copy_num",
            y="super_title",
            data=filt_copy_nums,
            jitter=True,
            color="0.0",
            alpha=0.5,
            ax=ax[0],
        )
        ax[0].set_xticks(ticks=range(0, int(max(filt_copy_nums["copy_num"]))+2))
        ax[0].set_xlim(xmin=-0.5, xmax=max(filt_copy_nums["copy_num"]) + 0.5)
        ax[0].set_ylabel("Super Population")
        ax[0].set_xlabel("Copy Number")
        sns.stripplot(
            x="copy_num",
            y="pop",
            data=filt_copy_nums,
            jitter=True,
            color="0.0",
            alpha=0.5,
            ax=ax[1],
        )
        ax[1].set_xticks(ticks=range(0, int(max(filt_copy_nums["copy_num"]))+2))
        ax[1].set_xlim(xmin=-0.5, xmax=max(filt_copy_nums["copy_num"]) + 0.5)
        ax[1].set_ylabel("Population")
        ax[1].set_xlabel("Copy Number")
        fig.subplots_adjust(wspace=0.4)
        ax[0].grid(True)
        ax[1].grid(True)

plt.savefig(args.output_name, format=args.file_type)
