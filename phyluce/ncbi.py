#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 09 June 2014 15:14 PDT (-0700)
"""

import os


def get_excludes(conf, sec):
    if conf.has_section(sec):
        excludes = [i[0] for i in conf.items(sec)]
    else:
        return []
    return excludes


def get_metadata(conf):
    return dict(conf.items("metadata"))


def get_vouchers(conf):
    if conf.has_section("vouchers"):
        return dict(conf.items("vouchers"))
    else:
        return None


def get_remaps(conf):
    if conf.has_section("remap"):
        return {i[0].replace(" ", "_"): i[1] for i in conf.items("remap")}
    else:
        return None


def get_species_name_from_file(infile, remap):
    sp = os.path.basename(infile).split(".")[0].replace("-", "_")
    return get_species_name(sp, remap)


def get_species_name(sp, remap):
    if remap is not None and sp in remap.keys():
        oldname = sp
        sp = remap[sp]
    else:
        oldname = sp
    species = sp.replace("_", " ").capitalize()
    partial = species.split(" ")[0].lower()[:3]
    return sp, species, partial, oldname


def get_node_name(read):
    # check for header match, if match get locus name for header
    nn_split = read.identifier.split("_")[:2]
    nn = "{}_{}".format(nn_split[0].strip(">").lower(), nn_split[1].lower())
    return nn


def get_new_identifier(species, uce, partial, counter, metadata, vouchers):
    title = "{0} ultra-conserved element locus {1}".format(species, uce)
    note = metadata["note"].format(uce)
    metadata = {
        "counter": counter,
        "partial": partial,
        "title": title,
        "organism": "{0}".format(species),
        "moltype": "{0}".format(metadata["moltype"]),
        "location": "{0}".format(metadata["location"]),
        "note": note,
        "specimen_voucher": "{0}".format(vouchers[species.lower()]),
    }
    # pdb.set_trace()
    new_identifier = "{counter}{partial} [organism={organism}] [moltype={moltype}] [location={location}] [note={note}] [specimen-voucher={specimen_voucher}] {title}".format(
        **metadata
    )
    return new_identifier
