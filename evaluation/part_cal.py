#!/usr/bin/env python
# coding=utf-8

def cal(static_dict, sigA, sigC, sigE, sigF, spacename, l):
    if spacename == "Total_aligned_bases":
        if sigA == True:
            static_dict["A"].Total_aligned_bases += l
        else:
            static_dict["B"].Total_aligned_bases += l

        if sigC:
            static_dict["C"].Total_aligned_bases += l
        else:
            static_dict["D"].Total_aligned_bases += l

        if sigE:
            static_dict["E"].Total_aligned_bases += l
        elif sigF:
            static_dict["F"].Total_aligned_bases += l
        else:
            static_dict["G"].Total_aligned_bases += l
    if spacename == "ExR80":
        if sigA:
            static_dict["A"].ExR80 += l
        else:
            static_dict["B"].ExR80 += l

        if sigC:
            static_dict["C"].ExR80 += l
        else:
            static_dict["D"].ExR80 += l

        if sigE:
            static_dict["E"].ExR80 += l
        elif sigF:
            static_dict["F"].ExR80 += l
        else:
            static_dict["G"].ExR80 += l

    if spacename == "ExA80":
        if sigA:
            static_dict["A"].ExA80 += l
        else:
            static_dict["B"].ExA80 += l

        if sigC:
            static_dict["C"].ExA80 += l
        else:
            static_dict["D"].ExA80 += l

        if sigE:
            static_dict["E"].ExA80 += l
        elif sigF:
            static_dict["F"].ExA80 += l
        else:
            static_dict["G"].ExA80 += l

    if spacename == "Hit80":
        if sigA:
            static_dict["A"].Hit80 += l
        else:
            static_dict["B"].Hit80 += l

        if sigC:
            static_dict["C"].Hit80 += l
        else:
            static_dict["D"].Hit80 += l

        if sigE:
            static_dict["E"].Hit80 += l
        elif sigF:
            static_dict["F"].Hit80 += l
        else:
            static_dict["G"].Hit80 += l


    if spacename == "ExR90":
        if sigA:
            static_dict["A"].ExR90 += l
        else:
            static_dict["B"].ExR90 += l

        if sigC:
            static_dict["C"].ExR90 += l
        else:
            static_dict["D"].ExR90 += l

        if sigE:
            static_dict["E"].ExR90 += l
        elif sigF:
            static_dict["F"].ExR90 += l
        else:
            static_dict["G"].ExR90 += l

    if spacename == "ExA90":
        if sigA:
            static_dict["A"].ExA90 += l
        else:
            static_dict["B"].ExA90 += l

        if sigC:
            static_dict["C"].ExA90 += l
        else:
            static_dict["D"].ExA90 += l

        if sigE:
            static_dict["E"].ExA90 += l
        elif sigF:
            static_dict["F"].ExA90 += l
        else:
            static_dict["G"].ExA90 += l

    if spacename == "ExR100":
        if sigA:
            static_dict["A"].ExR100 += l
        else:
            static_dict["B"].ExR100 += l

        if sigC:
            static_dict["C"].ExR100 += l
        else:
            static_dict["D"].ExR100 += l

        if sigE:
            static_dict["E"].ExR100 += l
        elif sigF:
            static_dict["F"].ExR100 += l
        else:
            static_dict["G"].ExR100 += l
    
    if spacename == "ExA100":
        if sigA:
            static_dict["A"].ExA100 += l
        else:
            static_dict["B"].ExA100 += l

        if sigC:
            static_dict["C"].ExA100 += l
        else:
            static_dict["D"].ExA100 += l

        if sigE:
            static_dict["E"].ExA100 += l
        elif sigF:
            static_dict["F"].ExA100 += l
        else:
            static_dict["G"].ExA100 += l
            
    if spacename == "Hit100":
        if sigA:
            static_dict["A"].Hit100 += l
        else:
            static_dict["B"].Hit100 += l

        if sigC:
            static_dict["C"].Hit100 += l
        else:
            static_dict["D"].Hit100 += l

        if sigE:
            static_dict["E"].Hit100 += l
        elif sigF:
            static_dict["F"].Hit100 += l
        else:
            static_dict["G"].Hit100 += l


    if spacename == "Total_aligned_exons":
        if sigA:
            static_dict["A"].Total_aligned_exons += l
        else:
            static_dict["B"].Total_aligned_exons += l

        if sigC:
            static_dict["C"].Total_aligned_exons += l
        else:
            static_dict["D"].Total_aligned_exons += l

        if sigE:
            static_dict["E"].Total_aligned_exons += l
        elif sigF:
            static_dict["F"].Total_aligned_exons += l
        else:
            static_dict["G"].Total_aligned_exons += l
    if spacename == "Total_aligned_reads":
        if sigA:
            static_dict["A"].Total_aligned_reads += l
        else:
            static_dict["B"].Total_aligned_reads += l

        if sigC:
            static_dict["C"].Total_aligned_reads += l
        else:
            static_dict["D"].Total_aligned_reads += l

        if sigE:
            static_dict["E"].Total_aligned_reads += l
        elif sigF:
            static_dict["F"].Total_aligned_reads += l
        else:
            static_dict["G"].Total_aligned_reads += l

