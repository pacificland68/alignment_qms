from django.shortcuts import render
from rest_framework.views import APIView
from rest_framework.response import Response
from rest_framework import status
from .models import alignment
import json

import Bio
import pandas as pd
import sys, getopt
import configparser
import os,re
from collections import defaultdict
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from urllib import request
from urllib.parse import urlparse
import requests
from django.http import FileResponse
import re

class AlignmentApiView(APIView):
    
    def post(self, req, *args, **kwargs):
        data = {

        }

        bin_path=os.path.abspath(os.path.dirname(sys.argv[0]))
        # clustal=bin_path+os.sep+"clustalw2.exe"
        clustal="clustalw2.exe"
        print(clustal)
        dir_cwd=os.getcwd()

        opts, args = getopt.getopt(sys.argv[1:],"hi:",["in="])
        if len(sys.argv) < 2:    
            print("python Sanger_Check_main.py -i <input.conf>")
            sys.exit(0)

        input_conf=os.path.join(os.path.dirname(__file__), 'input.ini')

        # for opt, arg in opts:
        #     print(arg, opt)
        #     if opt in ("-h","--help"):
        #         print("python Sanger_Check_main.py -i <input.conf>")
        #         sys.exit()
        #     elif opt in ("-i", "--in"):
        #         input_conf = arg

        conff=configparser.ConfigParser()
        conff.read(input_conf,encoding='UTF-8')
        sanger_dir=urlparse(conff.get("P1","sanger_dir"))
        outdir=conff.get("P1","outdir")
        # check_file=urlparse(conff.get("P1","check_file"))
        check_file=conff.get("P1","check_file")
        if os.path.exists(check_file):
            with open(check_file, 'rb') as file:
                response = FileResponse(file)
                print('response', response)

        if not os.path.exists(outdir):
            os.makedirs(outdir)

        fasta_dir=outdir+os.sep+"fasta"
        if not os.path.exists(fasta_dir):
            os.makedirs(fasta_dir)
        mut_dir=outdir+os.sep+"mutation_sites"
        #if not os.path.exists(mut_dir):
        #    os.makedirs(mut_dir)


        final_xls=outdir+os.sep+"final_result.xls"
        final_stat=outdir+os.sep+"Mutation_Ratio_stat.xls"
        #mutation_xls=mut_dir+os.sep+"Mutation_Site_result.xls"

        #########################################################################################################################
        ## multdict
        def nesteddict():
            return defaultdict(nesteddict)

        #########################################################################################################################
        ## Get diff num base of two sequence
        #def DiffBase(seq1,seq2,site_list):
        #    nsite=min(len(seq1),len(seq2))
        #    other_diff=0
        #    for n in range(nsite):
        #        if n in site_list:
        #            next
        #        else:
        #            base1=seq1[n]
        #            base2=seq2[n]
        #            if base1 != base2:
        #                other_diff=other_diff+1
        #    return other_diff



        #########################################################################################################################
        ## Get the Sanger start site
        def GetStartSite(file):
            align_dict=nesteddict()
            alignment = AlignIO.read(file, 'clustal')
            start=None
            for align in alignment:
                ID=align.id.split("-")[-1]
                if ID == "REF":
                    match=re.search(r"[A-Za-z]",str(align.seq))
                    start=match.span()[0]
                    ref_len=len(str(align.seq).replace("-",""))
                    for align2 in alignment:
                        TID=align2.id.split("-")[-1]
                        if TID != "REF":
                            seq=str(align2.seq)[start:].replace("-","")
                            if len(seq) > ref_len:
                                qseq=seq[:ref_len]
                            else:
                                qseq=seq
                            align_dict[align2.id]=qseq
            return align_dict

        ########################################################################################################################
        ## Get mutation Site
        #def GetMutationSite(seq1,seqid1,seq2,seqid2):



        #########################################################################################################################
        ## merge sanger seq to a fasta



        def GetSangerFile(path,odir):
            sangerf=nesteddict()
            print(path)
            for root,folders,files in os.walk(path):
                print('root', root)
                for file_name in files:
                    if file_name.endswith(".seq"):
                        file_p=root+os.sep+file_name
                        print(file_p)
                        Sanger_ID=re.findall("\(.*\)",file_name)
                        name_ID=Sanger_ID[0].strip("(").strip(")")
                        ID_len=name_ID.split("-")
                        if len(ID_len) > 3:
                            check_ID="-".join(name_ID.split("-")[1:-1])
                        else:
                            #len(name_ID) == 3:
                            check_ID="-".join(name_ID.split("-")[:-1])
                        #outfile=check_ID+".fa"
                        #outfile=odir+os.sep+outfile
                        #Sangerseq=open(file_p,'r')
                        Sangerseq=request.urlopen(file_p)
                        seq=Sangerseq.readline()
                        sangerf[check_ID][name_ID]=seq.decode()
                        print(str(check_ID)+"\t"+str(name_ID))
                        #with open(outfile,'a+') as F:
                            #F.write(">"+name_ID+"\n")
                            #F.write(seq+"\n")
                        #F.close()
            
            return sangerf


        print("Get The Sanger Files:\n")
        # SangerFile=GetSangerFile(sanger_dir,fasta_dir)
        print("\nChecking The Mutation Site ......\n")
        ##########################################################################################################################
        ## read check file
        # CF=open(check_file,'r')
        # CF=request.urlopen(check_file)
        CF=requests.get(check_file)
        CFR=CF.readlines()
        print('CFR[1]',CF.text)
        Final_P=open(final_xls,'w')
        #Detail_P=open(mutation_xls,'w')

        header=CFR[0].decode().strip("\n").split(",")
        final_header=CFR[0].decode().strip("\n").strip("\r").split(",")
        print('final_header', final_header)
        nhead=len(header)

        final_header.append("Sample\tMutation_Judge")

        Final_P.write("\t".join(final_header)+"\n")
        #Detail_P.write("\t".join(header)+"\n")



        ############################################################################################################################
        ## Check The Mutation Site in Sanger Sequence.
        stat_dict=nesteddict()
        Final_dict=nesteddict()
        Final_result=[]


        for i in range(1,len(CFR)):
            Frow=CFR[i].decode().strip("\n").split(",")
            check_id=Frow[0]+"-"+Frow[1]
            check_seq=Frow[3].upper()
            len_seq=len(check_seq)
            mutation_judge=nesteddict()
            nrow=nhead - len(Frow)
            #################### Get The Mutation Sites ################################################################################
            stat_dict[check_id]["True"]=0
            stat_dict[check_id]["False"]=0


            if nrow > 0:
                for n in range(nrow):
                    Frow.append(" ")

            Final_row=Frow
            #################### Get The Sanger Sequence ###############################################################################
            fafile=check_id+".fa"
            fa_trimmed=check_id+"_trimmed.fa"
            aligfile=check_id+".aln"
            #aalig_trimmed=check_id+"_trimmed.aln"
            FA_file=fasta_dir+os.sep+fafile
            trimmed=fasta_dir+os.sep+fa_trimmed
            Alig_file=fasta_dir+os.sep+aligfile
            #Alig_trimmed=fasta_dir+os.sep+alig_trimmed
            #print(str(check_id)+"\t"+str(check_seq))
            #Allseq=[]
            #Allseq.append(SeqRecord(Seq(check_seq), id=chec_id+"-Ref",annotations={"molecule_type": "DNA"}))
            if os.path.exists(FA_file):
                os.remove(FA_file)

            if check_id in SangerFile.keys():
                FA=open(FA_file,'a+')
                FA.write(">"+check_id+"-REF"+"\n")
                FA.write(check_seq+"\n")
                for k in SangerFile[check_id].keys():
                    FA.write(">"+k+"\n")
                    FA.write(SangerFile[check_id][k]+"\n")
                    #Allseq.append(SeqRecord(Seq(SangerFile[check_id][k]), id=k,annotations={"molecule_type": "DNA"}))

                FA.close()

            #################### do the alignment by clustalW ##############################################
            if os.path.exists(FA_file):
                os.chdir(fasta_dir)   
                clustal_cmd=ClustalwCommandline(clustal,infile=fafile,outfile=aligfile,align=True)
                assert os.path.isfile(clustal), "Clustal W not found"
                stdout, stderr = clustal_cmd()
                os.chdir(dir_cwd)
            
            ############################# Get The Start Site And Cut the Sanger Sequence ###################
            if os.path.exists(Alig_file):
                align_d=GetStartSite(Alig_file)
                
                if os.path.exists(trimmed):
                    os.remove(trimmed)
                FT=open(trimmed,'a+')
                FT.write(">"+check_id+"-REF"+"\n")
                FT.write(check_seq+"\n")
                
                for k in align_d.keys():
                    qseq=str(align_d[k])
                    FT.write(">"+k+"\n")
                    FT.write(qseq+"\n")
                    
                    if qseq == check_seq:
                        #Final_P.write("\t".join(Final_row)+"\t"+str(k)+"\t"+str("True")+"\n")
                        row="\t".join(Final_row)+"\t"+str(k)+"\t"+str("True")
                        Final_dict[check_id]["True"][k]=row
                        stat_dict[check_id]["True"]+=1
                    else:
                        #Final_P.write("\t".join(Final_row)+"\t"+str(k)+"\t"+str("False")+"\n")
                        row="\t".join(Final_row)+"\t"+str(k)+"\t"+str("False")
                        Final_dict[check_id]["False"][k]=row
                        stat_dict[check_id]["False"]+=1

                    

                FT.close()
            else:
                #Final_P.write("\t".join(Final_row)+"\t"+str(" ")+"\t"+str(" ")+"\n")
                row="\t".join(Final_row)+"\t"+str(" ")+"\t"+str(" ")
                Final_dict[check_id]["None"]["None"]=row



        ## end 
        ###################################################################################################################################################################
        for j in range(1,len(CFR)):
            Srow=CFR[j].decode().strip("\n").split(",")
            check_id=Srow[0]+"-"+Srow[1]
            tempResult=False
            if check_id in Final_dict.keys():
                if "True" in Final_dict[check_id].keys():
                    for sam in Final_dict[check_id]["True"].keys():
                        tempResult=True
                        # print(check_id,Final_dict[check_id]["True"][sam].split('\t')[16],Final_dict[check_id]["True"][sam].split('\t')[17])
                        Final_P.write(Final_dict[check_id]["True"][sam]+"\n")
                if "False" in Final_dict[check_id].keys():
                    for sam2 in Final_dict[check_id]["False"].keys():
                        Final_P.write(Final_dict[check_id]["False"][sam2]+"\n")
                if "None" in Final_dict[check_id].keys():
                    for sam3 in Final_dict[check_id]["None"].keys():
                        Final_P.write(Final_dict[check_id]["None"][sam3]+"\n")
            Final_result.append({
                "sample": check_id,
                "result": tempResult
            })
                    
        # print(Final_result)
        Final_P.close()
        ## end output 
        ###############################################################



        ST=open(final_stat,"w")
        ST.write("Sample\tTrue\tFalse\tTrue_Ratio%\n")
        for s in stat_dict.keys():
            total=stat_dict[s]["True"]+stat_dict[s]["False"]
            if total == 0:
                ST.write(str(s)+"\t0\t0\t0\n")
            else:
                true_ratio=stat_dict[s]["True"]/total
                true_ratio=round(true_ratio,4)*100
                ST.write(str(s)+"\t"+str(stat_dict[s]["True"])+"\t"+str(stat_dict[s]["False"])+"\t"+str(true_ratio)+"\n")

        ST.close()

        print("\nAll Sanger Files Have Been Checked.\n")

        return Response('here is alignment speaking', status=status.HTTP_200_OK )




class AlignmentNumericalApiView(APIView):

    def post(self, req, *args, **kwargs):
        # print('req.body', req.body)
        requestArray = json.loads(req.body)
        sampleParamResultList = []
        taskStepBatchId = requestArray["taskStepBatchId"]

        for item in requestArray["data"]:

            sampleId = item["sampleId"]
            sampleCode = item["sampleCode"]
            sampleName = item["sampleName"]
            
            for qcItem in item["itemList"]:
                # 单个检测项目下的所有检测参数
                for qcParam in qcItem["paramResultList"]:
                    result = 0
                    paramName = qcParam["paramName"]
                    greyRange = qcParam["greyArea"]
                    standarValue = qcParam["standardValue"]
                    # experimentalValue = qcParam["experimentalValue"]
                    experimentalValue = 1

                    standarChecker = RangeExpression(str(standarValue))
                    greyRangeChecker = RangeExpression(str(greyRange))
                
                    standarResult = standarChecker.in_range(experimentalValue)
                    greyResult = greyRangeChecker.in_range(experimentalValue)
                
                    # print(standarResult, greyResult)

                    if(not standarResult and not greyResult):
                        result = 2
                    elif greyResult:
                        result = 3
                    else: 
                        result = 1
                    
                    sampleParamResult = {
                        "actualValue": experimentalValue,
                        "greyArea":  greyRange,
                        "paramName": paramName,
                        "sampleId": sampleId,
                        "standardValue": standarValue,
                        "systemResult": result
                    }
                    sampleParamResultList.append(sampleParamResult)

        saveRequest = {
            "sampleParamResultList": sampleParamResultList,
            "taskStepBatchId": taskStepBatchId
        }

        # response = requests.post('https://wpis-dev.vpclub.cn/qms/inspectionRecord/saveParamResult', data=saveRequest)
        # if response.status_code == 200:
        #     saveResult = response.json()
        #     print('saveResult', saveResult)

        return Response(saveRequest, status=status.HTTP_200_OK)

class RangeExpression:
    def __init__(self, expression):
        self.ranges = self.parse_expression(expression)

    @staticmethod
    def parse_expression(expression):
        # 使用正则表达式匹配形如 [a,b] 和 [a,b) 的模式
        pattern = r'(\[|\()(\d+),(\d+)(\]|\))'
        matches = re.findall(pattern, expression)
        ranges = []
        
        for match in matches:
            start_inclusive = match[0] == '['
            end_inclusive = match[3] == ']'
            start_value = int(match[1])
            end_value = int(match[2])
            
            ranges.append((start_value, start_inclusive, end_value, end_inclusive))
        
        return ranges

    def in_range(self, num):
        for r in self.ranges:
            left_condition = (r[0] < num) if not r[1] else (r[0] <= num)
            right_condition = (r[2] > num) if not r[3] else (r[2] >= num)
            
            if left_condition and right_condition:
                return True
        return False


class Test(APIView):

    def post(self, req, *args, **kwargs):
        headers = {'Content-Type': 'application/json'}
        payload = {
            'id': 13
        }
        requestExperiment = {
            'batchTaskId': 13
        }
        qcItemResults = []

        # 根据批次id获取质控信息
        response = requests.post('https://wpis-dev.vpclub.cn/qms/inspectionRecord/getStepBatchInfo', data=payload)

        # 根据批次id获取样本实际值
        experimentResponse = requests.post('https://wpis-dev.vpclub.cn/pdm/projectBatchTaskSampleInspection/list', data=requestExperiment)

        if response.status_code == 200:
            qcItemResults = response.json()
            data = qcItemResults['data']

            requestParam = {
                "data": data,
                "taskStepBatchId": 13
            }

            # for sample in data:
            #     sampleId = sample["sampleId"]
            #     # 根据样本id获取实际值
            #     sampleCode = sample["sampleCode"]
            #     sampleName = sample["sampleName"]
            #     qcItemList = sample["itemList"]
            #     # 获取每个样本下的质控项目
            #     for qcItem in qcItemList:
            #         # 获取每个质控项目下的质控参数list
            #         qcParamList = qcItem["paramResultList"]
            #         for qcParam in qcParamList:
            #             # 计算每个质控参数的比对结果
            #             greyArea = qcParam["greyArea"]
            #             standardValue = qcItem["standardValue"]
            #             # 实际值
            #             experimentalValue = 1




        
        # payload = [
        #             {
        #                 "sampleId": 1,
        #                 "sampleCode": 2,
        #                 "sampleName": "material",
        #                 "paramResultList": [
        #                     {
        #                         "paramName": "",
        #                         "greyArea": "[1,2]",
        #                         "standarValue": "(2,3]",
        #                         "experimentalValue": 3
        #                     }
        #                 ]
        #             }
        #         ]
        response = requests.post('http://127.0.0.1:8000/alignment/api/number', json=requestParam, headers=headers)
        if response.status_code == 200:
            return Response(response.json(), status=status.HTTP_200_OK)
        else:
            return Response("", status=status.HTTP_200_OK)
        
        