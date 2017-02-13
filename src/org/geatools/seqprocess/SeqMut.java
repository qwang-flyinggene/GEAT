/*******************************************************************************
 *  ========================================================================
 *  GEATools : a free Genomic Event Analysis Tools for the Java(tm) platform
 *  ========================================================================
 *
 *  (C) Copyright 2016, by Qi Wang and Contributors.
 *
 *  This file is part of GEATools.
 *
 *  GEATools is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 *  [Oracle and Java are registered trademarks of Oracle and/or its affiliates. 
 *  Other names may be trademarks of their respective owners.]
 *
 *  (C) Copyright 2000-2014, by Original Author and Contributors.
 *
 *  Original Author: Qi Wang;
 *  Contributor(s): 
 *  Changes (from 1-Jan-2016)
 *  ---------------------------------------------------------------------------
 *******************************************************************************/
package org.geatools.seqprocess;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import org.geatools.data.structure.BaseMutStat;
import org.geatools.data.structure.BaseSNP;
import org.geatools.data.structure.BlastInfo;
import org.geatools.data.structure.SeqMutInfo;
import org.geatools.data.structure.SeqQual;
import org.geatools.data.structure.BaseMut;
import org.geatools.operation.FileOperate;

public class SeqMut {
    
	public static int colIdx_seqName=0;
	public static int colIdx_queryHitPos=4;
	public static int colIdx_refHitPos=5;
	public static int colIdx_strand=6;
	public static int colIdx_refBase=7;
	public static int colIdx_altBase=8;
	
	//int blastWordSize=7;
	//String blastTask="blastn-short";
	int blastWordSize=21;
	String blastTask="blastn";	
	String homeDir=".";
	String dataDir;
	List<String> tmpFiles;
	String tmpDir;	

	// Alignment Quality Standard(AQS) filter criteria	
	int minAlignLen=200;
	int maxMismatch=5;
	int maxGapNum=1;
	int maxQStart=50;
	int minQStart=1;
	int maxQEnd=500;
	int minQEnd=200;
	int maxSStart=50;
	int minSStart=1;
	int maxSEnd=500;
	int minSEnd=200;
	int totalSeqNum;
	int filteredSeqNum;
	
	//Neighbor Quality Standard(NQS) filter criteria
	int NQS_minS=30;
	int NQS_NN=5;
	int NQS_minNS=25;
	
	//Indel filter criteria
	int maxIndel=1;
	
	//Local Mutation Number Standard(LNS) filter criteria
	int localBaseNum=5;
	int localMaxMutNum=3;
	
	boolean doInterestMut=false;
	String [][]interestMut;
	
	float snpRate=0.2f;
	float baseErrRate=0.01f;
	
	List<BaseSNP> snpList;
	
	boolean skipAQS=false;
	boolean skipNQS=false;
	boolean skipIndel=false;
	boolean skipLMS=false;
	boolean skipBaseMutCheck=true;
	boolean doSNPCheck=false;
	boolean doBaseErrCheck=false;
	
	int avgRegion_start=1;
	int avgRegion_end=200;
	int refSeqStart;
	int refSeqEnd;
	
	boolean isReverse=false;
	
	public SeqMut(){
		
	}
	
	public void setHomeDir(String dir){
	    homeDir=dir;
	}
	  
	public void setDataDir(String dir){
		dataDir=dir;
	}
	  
	public void setTmpDir(String dir){
	    tmpDir=dir;
	}
	
	public void setAQS(int alignLen,int mismatch,int gapNum,boolean skip){
		minAlignLen=alignLen;
		maxMismatch=mismatch;
		maxGapNum=gapNum;
		skipAQS=skip;
	}
	
	public void setNQS(int NN,int minScore,int minNScore,boolean skip){
		NQS_NN=NN;
		NQS_minS=minScore;
		NQS_minNS=minNScore;
		skipNQS=skip;
	}
	
	public void setIndel(int indel,boolean skip){
	    maxIndel=indel;
		skipIndel=skip;
	}

	public void setLMS(int BN,int maxMN,boolean skip){
		localBaseNum=BN;
		localMaxMutNum=maxMN;
		skipLMS=skip;
	}
	
	public void setSNP(float rate,boolean doIt){
	    if(doIt) skipBaseMutCheck=false;
		snpRate=rate;
 	    doSNPCheck=doIt;
	}
	
	public void setBaseErr(float rate,boolean doIt){
		if(doIt) skipBaseMutCheck=false;
		baseErrRate=rate;
 	    doBaseErrCheck=doIt;
	}
	
	public void setTargetMut(String[][]mut){
		interestMut=mut;
		if(interestMut!=null && interestMut.length>0) doInterestMut=true;
		else doInterestMut=false;
	}
	
	public void setAvgRegion(int start, int end){
	  if(end>=start){
		avgRegion_start=start;
		avgRegion_end=end;
	  }else{
		avgRegion_start=end;
		avgRegion_end=start;
		System.out.println("Warning: your '-start' is large than '-end', the system automatically reverses them!");
	  }
	}
	
	public void setIsReversed(boolean status){
 	  isReverse=status;
	}
	
	public List<BaseSNP> getSNPList(){
		return snpList;
	}
	
	public SeqMutInfo getBLASTSeqMut(String querySeqFile, String refSeqFile,
			List<SeqQual>seqQualList){
	  
	    if(querySeqFile==null || refSeqFile==null) return null;
		  
		SeqMutInfo mutInfo=null;
		int res1=1;
		int res2=1;

		tmpFiles=new ArrayList<String>();
		
		if(doSNPCheck) snpList=new ArrayList<BaseSNP>();
		
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String blastOutFile=tmpDir+"/Mut.BlastOut."+timeStamp+".txt";
		String blastOutXMLFile=tmpDir+"/Mut.BlastOut."+timeStamp+".xml";
		tmpFiles.add(blastOutFile);
		tmpFiles.add(blastOutXMLFile);
		String blastCMD="blastn -task "+blastTask
		          +" -word_size="+blastWordSize
		          +" -max_target_seqs=1"
		          +" -query "+querySeqFile
		          +" -subject "+refSeqFile
		          +" -out "+blastOutFile
		          +" -outfmt 6";

		res1=SeqOperation.runBLAST(blastCMD);
		if(res1==0){
			System.out.println("Got BLAST txt!");
		}else{
			System.out.println("Error for BLAST txt out!!!");
			//System.exit(-1);
			return null;
		}
		
		blastCMD="blastn -task "+blastTask
		          +" -word_size="+blastWordSize
		          +" -max_target_seqs=1"
		          +" -query "+querySeqFile
		          +" -subject "+refSeqFile
		          +" -out "+blastOutXMLFile
		          +" -outfmt 5";
		res2=SeqOperation.runBLAST(blastCMD);
		if(res2==0){
			System.out.println("Got BLAST XML!");
		}else{
			System.out.println("Error for BLAST xml out!!!");
			//System.exit(-1);
			return null;
		}
		
		List<ArrayList<String>> blastOut=FileOperate.getMatrixFromFile(blastOutFile);
		totalSeqNum=blastOut.size();		
		List<BaseMut> baseMutList=getMutFromBLASTXML(blastOutXMLFile);
		System.out.println("Total sequences: "+totalSeqNum);	
		System.out.println("Total mutations: "+baseMutList.size());
		
		int usedSeqNum;
		List<BaseMut> filterBaseMutList = new ArrayList<BaseMut>();
		if(skipAQS){
		   System.out.println("skip AQS filter");
		   usedSeqNum=totalSeqNum;
		}else{		  
		   System.out.println("AQS filtering.........");
		   filterBaseMutList=filterBaseMut_AQS(baseMutList,blastOut);
		   blastOut=null;
		   usedSeqNum=filteredSeqNum;
		   System.out.println("Got sequences by AQS filter: "+filteredSeqNum);		   
		   System.out.println("Got mutations by AQS filter: "+filterBaseMutList.size());
		}
		
		if(skipNQS){
		   System.out.println("skip NQS filter");
		}else if(seqQualList!=null && seqQualList.size()>0){
		   System.out.println("NQS filtering.........");
		   filterBaseMutList=filterBaseMut_NQS(filterBaseMutList,seqQualList);
		   System.out.println("Got mutations by NQS filter: "+filterBaseMutList.size());
		}else{
		   System.out.println("Warning: No seq quality info, the system skip NQS");
		}
		
		if(skipIndel)
		   System.out.println("skip indel filter");
		else{
		   System.out.println("INDEL filtering.........");
		   filterBaseMutList=filterBaseMut_Indel(filterBaseMutList);
		   System.out.println("Got mutations by INDEL filter: "+filterBaseMutList.size());
		}
		
		if(skipLMS)
		   System.out.println("skip LMS filter");
		else{
		   System.out.println("LMS filtering.........");
		   filterBaseMutList=filterBaseMut_LMS(filterBaseMutList);
		   System.out.println("Got mutations by LMS filter: "+filterBaseMutList.size());
		}		
		
		mutInfo=getSeqMutInfo(filterBaseMutList,usedSeqNum,refSeqFile);	
		if(mutInfo!=null){
    	  mutInfo.expName=querySeqFile;
    	  mutInfo.seqName=refSeqFile;
          if(usedSeqNum==0){
            System.out.println("Warning: There is no any reads to pass filter criteria in "+querySeqFile);
          }else if(filterBaseMutList.size()==0){
            System.out.println("Warning: There is no any mutant base in "+querySeqFile);
          }else{

            System.out.println("Average BaseMutFrequency between "
                  +refSeqStart+"-"+refSeqEnd+" bp"
        		  +": "+mutInfo.avgMutRate + " for ["+querySeqFile+"]");
            if(doInterestMut){
            	System.out.println("Average BaseMutFrequency of interesting targets between "
                      +refSeqStart+"-"+refSeqEnd+" bp"
              		  +": "+mutInfo.avgTarMutRate + " for ["+querySeqFile+"]");
            	
            	System.out.println("Average BaseMutFrequency of others between "
                      +refSeqStart+"-"+refSeqEnd+" bp"
                	  +": "+mutInfo.avgOtherMutRate + " for ["+querySeqFile+"]");
            }
          }
		}else{
          System.out.println("Error: null  ["+refSeqFile+"]");
          return null;
        }
		filterBaseMutList=null;	
		baseMutList=null;
	
		for(String tmpFile: tmpFiles){
		    FileOperate.delFile(tmpFile);
		}

	    return mutInfo;
		
	}
	
	int getSeqNum(List<BaseMut> baseMutList){
		
		int seqNum=0;
		List<String> seqNameList=new ArrayList<String>();
		for(BaseMut mut: baseMutList){
			if(!seqNameList.contains(mut.seqName)) seqNameList.add(mut.seqName);
		}
		
		seqNum=seqNameList.size();
		seqNameList=null;
		
		return seqNum;
	}
	
	List <String> filterBLASTOut(List<ArrayList<String>> blastOut){		
	   
	   int alignLen=0;
	   int mismatchNum=0;
	   int gapNum=0;
	   int qStart=0;
	   //int qEnd=0;
	   //int sStart=0;
	   //int sEnd=0;
	   String seqName;		
	   List <String> seqNameList=new ArrayList <String>();
		
	   for(int i=0; i<blastOut.size();i++){
		 alignLen=Integer.parseInt(blastOut.get(i).get(BlastInfo.colAlignLen).trim());
		 mismatchNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colMismatchNum).trim());
		 gapNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colGapNum).trim());
		 qStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colQStart).trim());
		 //qEnd=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colQEnd).trim());
		 //sStart=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colSStart).trim());
		 //sEnd=Integer.parseInt(blastOut.get(i).get(blastInfoObj.colSEnd).trim());
		
		 if(alignLen>=minAlignLen && mismatchNum<=maxMismatch && gapNum<=maxGapNum 
					&& qStart<=maxQStart && qStart>=minQStart){
			 
			 seqName=blastOut.get(i).get(BlastInfo.colQName).trim();
			 seqNameList.add(seqName);
	
		}
	   }
		
	   return seqNameList;
		
	}
	
	List<BaseMut> filterBaseMut_AQS(List<BaseMut>baseMutList, 
			List<ArrayList<String>> blastOut){
        
		List<BaseMut>baseMutFilterList=new ArrayList<BaseMut>();
		List<String> filteredSeqName=filterBLASTOut(blastOut);
		filteredSeqNum=filteredSeqName.size();
		for(BaseMut baseMut:baseMutList){
			int idx=filteredSeqName.indexOf(baseMut.seqName);
			if(idx>=0) baseMutFilterList.add(baseMut);
		}
		
		return baseMutFilterList;
		
	}
	
	List<BaseMut> filterBaseMut_NQS(List<BaseMut>baseMutList,List<SeqQual>seqQualList){
        
		List<BaseMut>baseMutFilterList=new ArrayList<BaseMut>();
		if(seqQualList==null || seqQualList.size()==0) return baseMutList;
		if(baseMutList==null || baseMutList.size()==0) return baseMutFilterList;
		
		int queryPos;
		String seqName;
		int[] seqQual;
		int nStart=0;
		int nEnd=0;
		boolean isNeighborOk=false;		
		for(BaseMut mut:baseMutList){	
		   seqName=mut.seqName;
		   seqQual=getSeqQual(seqQualList,seqName);	
		   queryPos=mut.queryPos;
		   if(seqQual!=null && queryPos<seqQual.length){		
			  if(seqQual[queryPos-1]>=NQS_minS){				
				nStart=Math.max(0, queryPos-1-NQS_NN);
				nEnd=Math.min(seqQual.length-1, queryPos-1+NQS_NN);
				isNeighborOk=true;
				for(int i=nStart;i<nEnd;i++){
				  if(i!=(queryPos-1) && seqQual[i]<NQS_minNS){
					  isNeighborOk=false;
					  break;
				  }
				}
				if(isNeighborOk) baseMutFilterList.add(mut);
			  }
		   }			
		}
		
		return baseMutFilterList;
		
	}
	
	int[] getSeqQual(List<SeqQual>seqQualList, String seqName){
		
		int[] seqQual = null;
		for(SeqQual qual:seqQualList){	
			if(seqName.equalsIgnoreCase(qual.seqName)){
				seqQual=qual.seqQual;			
				break;
			}
		}
		
		return seqQual;
	}
	
	List<BaseMut> filterBaseMut_Indel(List<BaseMut>baseMutList){

		List<BaseMut>baseMutFilterList=new ArrayList<BaseMut>();
		int indel=0;
		for(BaseMut baseMut:baseMutList){
		  indel=0;
		  
		  if(baseMut.blastQueryBase.indexOf("-")>=0) 
			  indel=baseMut.blastQueryBase.length();		    
		  else if(baseMut.blastHitBase.indexOf("-")>=0) 
			  indel=baseMut.blastHitBase.length();
		  
		  if(indel<=maxIndel) baseMutFilterList.add(baseMut);
		}
		
		return baseMutFilterList;
		
	}
		
	
	List<BaseMut> filterBaseMut_LMS(List<BaseMut>baseMutList){

		List<BaseMut>baseMutFilterList=new ArrayList<BaseMut>();
		
		if(baseMutList==null || baseMutList.size()==0) return baseMutFilterList;
		
		BaseMut mut=baseMutList.get(0);
		String seqName=mut.seqName;		
		int queryPos=mut.queryPos;	
		int mutNum=getBaseMutNum(mut);	
		
		BaseMut mut1;
		String seqName1;
		int queryPos1;
		int mutNum1=0;

		int i=1;
		while(i<baseMutList.size()){
		
		  mut1=baseMutList.get(i);
		  seqName1=mut1.seqName;
		  queryPos1=mut1.queryPos;
		  mutNum1=getBaseMutNum(mut1);
		  if(seqName.equalsIgnoreCase(seqName1)){		
			 if((queryPos1-queryPos)>localBaseNum){			 	
				if(mutNum<=localMaxMutNum) baseMutFilterList.add(mut);	
				mutNum=mutNum1;				
			    seqName=mut1.seqName;		
				queryPos=mut1.queryPos;
				mut=mut1;
			 }else{
			    mutNum=mutNum+mutNum1;			   
			 }

			 i++;		
		  }else{
			 if(mutNum<=localMaxMutNum) baseMutFilterList.add(mut);	
			 mut=mut1;
			 seqName=mut.seqName;		
			 queryPos=mut.queryPos;	
			 i++;
		  }

		}
		
		return baseMutFilterList;
		
	}
	
    int getBaseMutNum(BaseMut mutObj){
		int mutNum=0;
		String refBase;
		String altBase;
		int del=mutObj.blastQueryBase.indexOf("-");
		int insert=mutObj.blastHitBase.indexOf("-");
		if(del<0 && insert<0){ 
			mutNum=0;
			refBase=mutObj.refBase;
			altBase=mutObj.altBase;
		    for(int j=0;j<altBase.length();j++){
			   if(refBase.charAt(j)!=altBase.charAt(j)){
				   mutNum=mutNum+1;		          
			   }
			}		
		}
		
		return mutNum;
	}
    
    int getBaseMutNum(BaseMut mutObj,String [][]mutBasePairs){
		int mutNum=0;
		String refBase;
		String altBase;
		String rBasei;
		String aBasei;
		int del=mutObj.blastQueryBase.indexOf("-");
		int insert=mutObj.blastHitBase.indexOf("-");
		if(del<0 && insert<0){ 
		   mutNum=0;
		   refBase=mutObj.refBase;
		   altBase=mutObj.altBase;
		   if(mutBasePairs!=null && mutBasePairs.length>0){
		     for(int i=0;i<altBase.length();i++){
		       rBasei=refBase.substring(i,(i+1));
		       aBasei=altBase.substring(i,(i+1));
		       for(int m=0;m<mutBasePairs.length;m++){
			     if(!rBasei.equalsIgnoreCase(aBasei) 
			    	  && rBasei.equalsIgnoreCase(mutBasePairs[m][0]) 
			    	  && aBasei.equalsIgnoreCase(mutBasePairs[m][1])){
			    	 
				    mutNum=mutNum+1;
				    break;				   
			     }
		       }
			 }
		   }else{
			 for(int i=0;i<altBase.length();i++){
			   if(refBase.charAt(i)!=altBase.charAt(i)){
				   mutNum=mutNum+1;		          
			   }
			 }		 
		   }
		}
		
		return mutNum;
	}
	
	BaseMut getBaseMutByQueryPos(List<BaseMut> baseMutList, int queryPos){
	   BaseMut baseMut = null;
	   for(BaseMut mut:baseMutList){
		   if(queryPos==mut.queryPos){
			   baseMut=mut;
			   break;
		   }
	   }
	   
	   return baseMut;
	}
	
	List<BaseMut> getMutFromBLASTXML(String blastOutXMLFile){
		
		String []args=new String[3];
		args[0]="-o";
		args[1]=blastOutXMLFile+".mut";	
		args[2]=blastOutXMLFile;
		//new BlastNToSnp().instanceMainWithExit(args);	
		BlastNToSnp blastn2mut=new BlastNToSnp();
		int res=blastn2mut.doWork(args);	
		blastn2mut=null;
		
		int rowSize;
		String seqName;
		List<BaseMut> baseMutList=new ArrayList<BaseMut>();
		BaseMut baseMut;	
		if(res==0){
		   List<ArrayList<String>> blast2mut=FileOperate.getMatrixFromFile(args[1]); 		
		   for(int i=1;i<blast2mut.size();i++){
			  baseMut=new BaseMut();
			  rowSize=blast2mut.get(i).size();
			  seqName=blast2mut.get(i).get(colIdx_seqName).trim();	
			  //System.out.println(seqName);
		      if(seqName.indexOf(" ")>0)
			    baseMut.seqName=seqName.substring(0,seqName.indexOf(" "));
		      else
		    	baseMut.seqName=seqName;
			  
		      baseMut.queryPos=Integer.parseInt(blast2mut.get(i).get(colIdx_queryHitPos).trim());
			  baseMut.refPos=Integer.parseInt(blast2mut.get(i).get(colIdx_refHitPos).trim());
		      baseMut.altBase=blast2mut.get(i).get(colIdx_altBase).trim();
			  baseMut.refBase=blast2mut.get(i).get(colIdx_refBase).trim();
			  baseMut.strand=blast2mut.get(i).get(colIdx_strand).trim();
			  baseMut.blastQueryBase=blast2mut.get(i).get(rowSize-1-1).trim();
			  baseMut.blastHitBase=blast2mut.get(i).get(rowSize-1-2).trim();
			  if(baseMut.blastQueryBase.indexOf("-")>=0) 
				  baseMut.altBase=baseMut.blastQueryBase;
			  if(baseMut.blastHitBase.indexOf("-")>=0) 
				  baseMut.refBase=baseMut.blastHitBase;
			  
				  
			  baseMutList.add(baseMut);
		   }
		}else{
			System.exit(res);
		}
		
		return baseMutList;
		
	}
	
	public SeqMutInfo getSeqMutInfo(List<BaseMut> baseMutList, int readNum,
			String refSeqFastaFile){
		
		if(!SeqOperation.isFASTASeq(refSeqFastaFile)) return null;
		
		int refHitPos=0;	
	
		SeqMutInfo mutInfo=new SeqMutInfo();
		mutInfo.readsNum=readNum;		
		mutInfo.seq=SeqOperation.getFASTASeqObj(refSeqFastaFile).get(0).seq;
		mutInfo.baseMutStat=new BaseMutStat[mutInfo.seq.length()];
		for(int i=0;i<mutInfo.seq.length();i++){
			mutInfo.baseMutStat[i]=new BaseMutStat();
			mutInfo.baseMutStat[i].mutCount=0;
			mutInfo.baseMutStat[i].mutRate=0;
			mutInfo.baseMutStat[i].tarMutCount=0;
			mutInfo.baseMutStat[i].tarMutRate=0;
			mutInfo.baseMutStat[i].otherMutCount=0;
			mutInfo.baseMutStat[i].otherMutRate=0;
			mutInfo.baseMutStat[i].a=0;
			mutInfo.baseMutStat[i].g=0;
			mutInfo.baseMutStat[i].c=0;
			mutInfo.baseMutStat[i].t=0;
			mutInfo.baseMutStat[i].del=0;
			mutInfo.baseMutStat[i].insert=0;		
		}
		
		if(baseMutList==null || baseMutList.size()==0) return mutInfo;
		
		String altBase;
		String refBase;
		String aBasei;
		String rBasei;
		for(BaseMut baseMut:baseMutList){		
		  altBase=baseMut.altBase;
		  refBase=baseMut.refBase;
		  refHitPos=baseMut.refPos;
		  if(baseMut.blastHitBase.indexOf("-")<0){ //exclude inserted bases in seq			   		    
			int alignRefStartIdx=refHitPos-1;
			if(baseMut.strand.trim().equals("-")) alignRefStartIdx=refHitPos-altBase.length();
			for(int i=0;i<altBase.length();i++){
				rBasei=String.valueOf(refBase.charAt(i));
				aBasei=String.valueOf(altBase.charAt(i));
				if(!rBasei.equalsIgnoreCase(aBasei)){				        
				     if(aBasei.equalsIgnoreCase("A"))
				          mutInfo.baseMutStat[alignRefStartIdx+i].a++;
				     if(aBasei.equalsIgnoreCase("G"))
					      mutInfo.baseMutStat[alignRefStartIdx+i].g++;
				     if(aBasei.equalsIgnoreCase("C"))
					      mutInfo.baseMutStat[alignRefStartIdx+i].c++;
				     if(aBasei.equalsIgnoreCase("T"))
					      mutInfo.baseMutStat[alignRefStartIdx+i].t++;
				     if(aBasei.equalsIgnoreCase("-"))
						  mutInfo.baseMutStat[alignRefStartIdx+i].del++;
				     
				     mutInfo.baseMutStat[alignRefStartIdx+i].mutCount++;				 
				}
			 }		    
		  }else{
			//count inserted mutant bases into refBase where start to insert	
			 mutInfo.baseMutStat[refHitPos-1].insert=altBase.length();
			 mutInfo.baseMutStat[refHitPos-1].mutCount
					=mutInfo.baseMutStat[refHitPos-1].mutCount+altBase.length(); 	
		  }
		}
		refSeqStart=avgRegion_start;
		refSeqEnd=avgRegion_end;
		if(isReverse){
			refSeqStart=mutInfo.seq.length()-avgRegion_end+1;
			refSeqEnd=mutInfo.seq.length()-avgRegion_start+1;	
		}
		double sumMutRate=0;
		int n=0;
		if(skipBaseMutCheck){
		  n=(refSeqEnd-refSeqStart+1);
		  for(int i=0;i<mutInfo.seq.length();i++){			
			mutInfo.baseMutStat[i].mutRate=(1.0d*mutInfo.baseMutStat[i].mutCount)/mutInfo.readsNum;				
			if(i>=refSeqStart && i<=refSeqEnd){
				sumMutRate=sumMutRate+mutInfo.baseMutStat[i].mutRate;		
			}
		  }
		}else{			
		  snpList=new ArrayList<BaseSNP>();		
		  for(int i=0;i<mutInfo.seq.length();i++){			
			 mutInfo.baseMutStat[i].mutRate=(1.0d*mutInfo.baseMutStat[i].mutCount)/mutInfo.readsNum;	
			 if(doSNPCheck &&((1.0d*mutInfo.baseMutStat[i].g/mutInfo.readsNum)>=snpRate
				 || (1.0d*mutInfo.baseMutStat[i].c/mutInfo.readsNum)>=snpRate
				 || (1.0d*mutInfo.baseMutStat[i].t/mutInfo.readsNum)>=snpRate
				 || (1.0d*mutInfo.baseMutStat[i].a/mutInfo.readsNum)>=snpRate)){
				 
				  BaseSNP snp=new BaseSNP();	
				  BaseMut base=getBaseMut(baseMutList,(i+1));
				  if(base!=null){										
					snp.refBase=base.refBase;
					snp.refPos=base.refPos;					
					base=null;
				  }
				  snp.altBaseStat=mutInfo.baseMutStat[i];
				  snp.snpRate=mutInfo.baseMutStat[i].mutRate;
				  snpList.add(snp);
				  snp=null;
				  
				  mutInfo.baseMutStat[i].mutRate=-1;
			 }
			 
			 if(doBaseErrCheck &&(mutInfo.baseMutStat[i].mutRate>=baseErrRate)){
				  mutInfo.baseMutStat[i].mutRate=-2;
			 }
			 
			 if(i>=refSeqStart && i<=refSeqEnd && mutInfo.baseMutStat[i].mutRate>=0){
				 sumMutRate=sumMutRate+mutInfo.baseMutStat[i].mutRate;	
				 n=n+1;
			 }
		  }// for i
		}

	    mutInfo.sumMutRate=sumMutRate;
		mutInfo.baseNum=n;
		mutInfo.avgMutRate=sumMutRate/n;
		
		if(doInterestMut) setTargetMutInfo(mutInfo);
		
		return mutInfo;
	}
	
	BaseMut getBaseMut(List<BaseMut>baseMutList, int refPos){
		
		BaseMut base = null;
		for(BaseMut baseMut:baseMutList){
		   if(refPos==baseMut.refPos){	
				base=baseMut;
				break;					  	
		   }
	    }
		
		return base;
	}
	
	void setTargetMutInfo(SeqMutInfo mutInfo){
		

		if(interestMut!=null && interestMut.length>0){
		   String refBasei;	
		   double sumTarMutRate=0.0d;
		   double sumOtherMutRate=0.0d;		 
		   int targetBaseNum=0;
		   int otherBaseNum=0;
		   int tarMutCount=0;
		   boolean isValidTarBase=false;
		   boolean isValidBase=false;
		   String[] itemSplited;
		   boolean isGUsed=false;
		   boolean isCUsed=false;
		   boolean isTUsed=false;
		   boolean isAUsed=false;
		   boolean isDelUsed=false;
		   boolean isNUsed=false;
		   for(int i=0;i<mutInfo.seq.length();i++){
			  isValidBase=false;
			  isValidTarBase=false;	
			  tarMutCount=0;
			  mutInfo.baseMutStat[i].tarMutCount=0;			  
			  /*
			  mutInfo.baseMutStat[i].otherMutCount=   
			             +mutInfo.baseMutStat[i].g
					     +mutInfo.baseMutStat[i].c
                         +mutInfo.baseMutStat[i].a
                         +mutInfo.baseMutStat[i].t
                         +mutInfo.baseMutStat[i].del;
			  */
			  if(i>=refSeqStart && i<=refSeqEnd) isValidBase=true;	
			  refBasei=String.valueOf(mutInfo.seq.charAt(i));
			  for(int m=0;m<interestMut.length;m++){
			     if(refBasei.equalsIgnoreCase(interestMut[m][0])){
			
			    	if(i>=refSeqStart && i<=refSeqEnd) isValidTarBase=true;			    		 
			    	 /*
			    	if(interestMut[m][1].equalsIgnoreCase("G")){			    		   
					   	 mutInfo.baseMutStat[i].tarMutCount=mutInfo.baseMutStat[i].g;					   
					   	 break;
				    }else if(interestMut[m][1].equalsIgnoreCase("C")){			    		   
						 mutInfo.baseMutStat[i].tarMutCount=mutInfo.baseMutStat[i].c;					   	 
					     break;
				    }else if(interestMut[m][1].equalsIgnoreCase("A")){			    		   
					   	 mutInfo.baseMutStat[i].tarMutCount=mutInfo.baseMutStat[i].a;					
					   	 break;
			    	}else if(interestMut[m][1].equalsIgnoreCase("T")){			    		   
						 mutInfo.baseMutStat[i].tarMutCount=mutInfo.baseMutStat[i].t;					   
					     break;
			    	}else if(interestMut[m][1].equalsIgnoreCase("-")){			    		   
						 mutInfo.baseMutStat[i].tarMutCount=mutInfo.baseMutStat[i].del;					   	
						 break;
				    }else if(interestMut[m][1].equalsIgnoreCase("N")){			    		   
						 mutInfo.baseMutStat[i].tarMutCount=mutInfo.baseMutStat[i].mutCount;
						 break;
					}else{			    	  
			    	    break;
					}
					*/
			    	
			 	    isGUsed=false;
				    isCUsed=false;
				    isTUsed=false;
				    isAUsed=false;
				    isDelUsed=false;
				    isNUsed=false;
				    mutInfo.baseMutStat[i].tarMutCount=0;					
			    	itemSplited=interestMut[m][1].split("|");
			    	for(String mut:itemSplited){
			    	   tarMutCount=0;
			    	   if(mut.equalsIgnoreCase("G") && !isGUsed){			    		   
			    		  tarMutCount=mutInfo.baseMutStat[i].g;
			    		  isGUsed=true;
					   }else if(mut.equalsIgnoreCase("C") && !isCUsed){			    		   
						  tarMutCount=mutInfo.baseMutStat[i].c;
						  isCUsed=true;
					   }else if(mut.equalsIgnoreCase("A") && !isAUsed){			    		   
						  tarMutCount=mutInfo.baseMutStat[i].a;	
						  isAUsed=true;
				       }else if(mut.equalsIgnoreCase("T") && !isTUsed){			    		   
				    	  tarMutCount=mutInfo.baseMutStat[i].t;	
				    	  isTUsed=true;
				       }else if(mut.equalsIgnoreCase("-") && !isDelUsed){			    		   
				    	  tarMutCount=mutInfo.baseMutStat[i].del;
				    	  isDelUsed=true;
					   }else if(mut.equalsIgnoreCase("N") && !isNUsed){			    		   
						  tarMutCount=mutInfo.baseMutStat[i].mutCount;
						  isNUsed=true;
					   }			    	  
			    	   mutInfo.baseMutStat[i].tarMutCount
			    	       =mutInfo.baseMutStat[i].tarMutCount+tarMutCount;
			    	}
			    	break;
				 }
			  }
			 	
			  mutInfo.baseMutStat[i].tarMutRate
					 =(1.0d*mutInfo.baseMutStat[i].tarMutCount)/mutInfo.readsNum;
			  
			  mutInfo.baseMutStat[i].otherMutCount
			         =mutInfo.baseMutStat[i].mutCount-mutInfo.baseMutStat[i].tarMutCount;
			  mutInfo.baseMutStat[i].otherMutRate
				     =(1.0d*mutInfo.baseMutStat[i].otherMutCount)/mutInfo.readsNum;
			
			  if(!skipBaseMutCheck && doBaseErrCheck){			    
			    if(mutInfo.baseMutStat[i].mutRate<0 
			    		|| mutInfo.baseMutStat[i].mutRate>=baseErrRate){
					
					isValidTarBase=false;
					isValidBase=false;

					if(mutInfo.baseMutStat[i].tarMutRate>=baseErrRate)
					   mutInfo.baseMutStat[i].tarMutRate=-2;
					if(mutInfo.baseMutStat[i].otherMutRate>=baseErrRate)
					   mutInfo.baseMutStat[i].otherMutRate=-2;
					if(mutInfo.baseMutStat[i].mutRate>=baseErrRate)
					   mutInfo.baseMutStat[i].mutRate=-2;
			    }
			  }
			  
			  if(isValidTarBase){
				  sumTarMutRate=sumTarMutRate+mutInfo.baseMutStat[i].tarMutRate;
				  targetBaseNum++;
			  }
			  if(isValidBase){
				  sumOtherMutRate=sumOtherMutRate+mutInfo.baseMutStat[i].otherMutRate;
				  otherBaseNum++;
			  }			 
		   }
		   mutInfo.sumTarMutRate=sumTarMutRate;
		   mutInfo.tarBaseNum=targetBaseNum;
		   mutInfo.avgTarMutRate=sumTarMutRate/targetBaseNum;
		   
		   mutInfo.sumOtherMutRate=sumOtherMutRate;
		   mutInfo.otherBaseNum=otherBaseNum;
		   mutInfo.avgOtherMutRate=sumOtherMutRate/otherBaseNum;
	   }
	}
	
	
	public void saveMutInfo(SeqMutInfo mutInfo, String outDir, String outTag){
		
		if(mutInfo==null) return;
		
		List<ArrayList<String>> out=new ArrayList<ArrayList<String>>();
		ArrayList<String> perRow;
		
		perRow=new ArrayList<String>();
		perRow.add("Base");
		perRow.add("Rate_Mut:All");
		perRow.add("Count_Mut:All");
		if(doInterestMut){
			String mutTag="";
			for(int m=0;m<interestMut.length;m++){
			  mutTag=mutTag+interestMut[m][0]+">"+interestMut[m][1]+"|";
			}
			mutTag=mutTag.substring(0,mutTag.lastIndexOf("|"));
			perRow.add("Rate_Mut:"+mutTag);
			perRow.add("Count_Mut:"+mutTag);
			perRow.add("Rate_Mut:Other");
			perRow.add("Count_Mut:Other");
		}
		perRow.add("G");
		perRow.add("C");
		perRow.add("A");
		perRow.add("T");
		perRow.add("Deletion");
		perRow.add("Insertion");
		out.add(perRow);
		for(int i=0;i<mutInfo.seq.length();i++){
			perRow=new ArrayList<String>();
			perRow.add(String.valueOf(mutInfo.seq.charAt(i)));
			if(mutInfo.baseMutStat[i].mutRate>=0)				
			  perRow.add(Double.toString(mutInfo.baseMutStat[i].mutRate));
			else
			  perRow.add("NA");
			perRow.add(Integer.toString(mutInfo.baseMutStat[i].mutCount));
			if(doInterestMut){
			  if(mutInfo.baseMutStat[i].tarMutRate>=0)
				perRow.add(Double.toString(mutInfo.baseMutStat[i].tarMutRate));	
			  else
				perRow.add("NA");
			  
			  perRow.add(Integer.toString(mutInfo.baseMutStat[i].tarMutCount));
			  
			  if(mutInfo.baseMutStat[i].otherMutRate>=0)
				perRow.add(Double.toString(mutInfo.baseMutStat[i].otherMutRate));	
			  else
				perRow.add("NA");
			  
			  perRow.add(Integer.toString(mutInfo.baseMutStat[i].otherMutCount));
			}
			perRow.add(Integer.toString(mutInfo.baseMutStat[i].g));
			perRow.add(Integer.toString(mutInfo.baseMutStat[i].c));
			perRow.add(Integer.toString(mutInfo.baseMutStat[i].a));
			perRow.add(Integer.toString(mutInfo.baseMutStat[i].t));
			perRow.add(Integer.toString(mutInfo.baseMutStat[i].del));
			perRow.add(Integer.toString(mutInfo.baseMutStat[i].insert));
			out.add(perRow);
			perRow=null;
		}
		
		String outFile;
		if(outDir!=null && new File(outDir).exists())
		  outFile=outDir+"/BaseMutInfo-"+outTag+".txt";
		else
		  outFile="BaseMutInfo-"+outTag+".txt";
		
		FileOperate.saveMatrixList(out, outFile);
		out=null;
	}
	
    public void saveBaseMutRate(List<SeqMutInfo>mutList,String outDir,String outTag){
		
    	if(mutList==null || mutList.size()<1) return;
		
		List<ArrayList<String>> out=new ArrayList<ArrayList<String>>();
		ArrayList<String> perRow;
		
		//Save all base mut frequency for each single base as profile
		out=new ArrayList<ArrayList<String>>();
		perRow=new ArrayList<String>();
		perRow.add("Pos");
		perRow.add("Base");	
	    for(SeqMutInfo mut:mutList){
	      perRow.add(mut.expName.substring(mut.expName.lastIndexOf("/")+1,mut.expName.length()));		
		}
	    out.add(perRow);
		perRow=null;
		  
		SeqMutInfo mut0=mutList.get(0);		
	    for(int i=0;i<mut0.seq.length();i++){
		  perRow=new ArrayList<String>();
		  perRow.add(Integer.toString(i+1));
		  perRow.add(String.valueOf(mut0.seq.charAt(i)));	
	      for(SeqMutInfo mut:mutList){
	    	if(mut.baseMutStat[i].mutRate>=0)	
			  perRow.add(Double.toString(mut.baseMutStat[i].mutRate));
	    	else
	    	  perRow.add("NA");	
		  }
		  out.add(perRow);
		  perRow=null;
		}
		String outFile;	
		if(outDir!=null && new File(outDir).exists())
		  outFile=outDir+"/BaseMutFrequency_"+outTag+".txt";
		else
		  outFile="BaseMutFrequency_"+outTag+".txt";
		
		FileOperate.saveMatrixList(out, outFile);
		out=null;
		
		//Save interesting target mut rate for each single base as profile
		if(doInterestMut){
			//saver intersting target base mut frequency
			out=new ArrayList<ArrayList<String>>();
			perRow=new ArrayList<String>();
			perRow.add("Pos");
			perRow.add("Base");	
		    for(SeqMutInfo mut:mutList){		    
			  perRow.add(mut.expName.substring(mut.expName.lastIndexOf("/"),mut.expName.length()));	
			}
		    out.add(perRow);
			perRow=null;
			  
			mut0=mutList.get(0);		
		    for(int i=0;i<mut0.seq.length();i++){
			  perRow=new ArrayList<String>();
			  perRow.add(Integer.toString(i+1));
			  perRow.add(String.valueOf(mut0.seq.charAt(i)));	
		      for(SeqMutInfo mut:mutList){
		    	if(mut.baseMutStat[i].tarMutRate>=0)	
				  perRow.add(Double.toString(mut.baseMutStat[i].tarMutRate));
		    	else
		    	  perRow.add("NA");
			  }
			  out.add(perRow);
			  perRow=null;
			}
					
			if(outDir!=null && new File(outDir).exists())
			  outFile=outDir+"/TargetBaseMutFrequency_"+outTag+".txt";
			else
			  outFile="TargetBaseMutFrequency_"+outTag+".txt";
			
			FileOperate.saveMatrixList(out, outFile);
			out=null;
			
			// save other base mut frequency
			out=new ArrayList<ArrayList<String>>();
			perRow=new ArrayList<String>();
			perRow.add("Pos");
			perRow.add("Base");	
		    for(SeqMutInfo mut:mutList){		    
			  perRow.add(mut.expName.substring(mut.expName.lastIndexOf("/"),mut.expName.length()));	
			}
		    out.add(perRow);
			perRow=null;
			  
			mut0=mutList.get(0);		
		    for(int i=0;i<mut0.seq.length();i++){
			  perRow=new ArrayList<String>();
			  perRow.add(Integer.toString(i+1));
			  perRow.add(String.valueOf(mut0.seq.charAt(i)));	
		      for(SeqMutInfo mut:mutList){
		    	if(mut.baseMutStat[i].otherMutRate>=0)	
				  perRow.add(Double.toString(mut.baseMutStat[i].otherMutRate));
		    	else
		    	  perRow.add("NA");
			  }
			  out.add(perRow);
			  perRow=null;
			}
					
			if(outDir!=null && new File(outDir).exists())
			  outFile=outDir+"/OtherBaseMutFrequency_"+outTag+".txt";
			else
			  outFile="OtherBaseMutFrequency_"+outTag+".txt";
			
			FileOperate.saveMatrixList(out, outFile);
			out=null;
		}
	}
    
    public void saveAvgBaseMutRate(List<SeqMutInfo>mutList,String outDir, String outTag){
		
    	if(mutList==null || mutList.size()<1) return;
    	
		List<ArrayList<String>> out=new ArrayList<ArrayList<String>>();
		ArrayList<String> perRow;
		
		//Save average base mut rate as profile
		perRow=new ArrayList<String>();
		perRow.add("Item");	
	    for(SeqMutInfo mut:mutList){	
	      perRow.add(mut.expName.substring(mut.expName.lastIndexOf("/")+1,mut.expName.length()));		
		}
	    out.add(perRow);
		perRow=null;
		
		//save reads count for each sample
	    perRow=new ArrayList<String>();
	    perRow.add("Reads Count");
	    for(SeqMutInfo mut:mutList){
	      perRow.add(Integer.toString(mut.readsNum));
		}
	    out.add(perRow);
		perRow=null;
		
		//save Average Base Mut
		if(doInterestMut){
		  String mutTag="";
		  for(int m=0;m<interestMut.length;m++){
			  mutTag=mutTag+interestMut[m][0]+">"+interestMut[m][1]+"|";
		  }
		  mutTag=mutTag.substring(0,mutTag.lastIndexOf("|"));
		  perRow=new ArrayList<String>();
		  perRow.add("Average BMF:"+mutTag);
	      for(SeqMutInfo mut:mutList){
	    	 perRow.add(Double.toString(mut.avgTarMutRate));
		  }
	      out.add(perRow);
		  perRow=null;
			
	      perRow=new ArrayList<String>();
	      perRow.add("Average BMF:Other");
	      for(SeqMutInfo mut:mutList){
	    	 perRow.add(Double.toString(mut.avgOtherMutRate));
		  }
	      out.add(perRow);
		  perRow=null;
		}		
		
	    perRow=new ArrayList<String>();
	    perRow.add("Average BMF:All");
	    for(SeqMutInfo mut:mutList){
	      perRow.add(Double.toString(mut.avgMutRate));
		}
	    out.add(perRow);
		perRow=null;
		
		String outFile="";
		if(outTag==null) outTag="profile";
				
		if(outDir!=null && new File(outDir).exists())
		  outFile=outDir+"/AverageBaseMutFrequency-"+outTag+".txt";
		else
		  outFile="AverageBaseMutFrequency-"+outTag+".txt";
		
		FileOperate.saveMatrixList(out, outFile);
		out=null;
	}
    
    public void saveAvgBaseMutRate(List<SeqMutInfo> mutList,List<SeqMutInfo> mutList2,
    		String outDir, String outTag){
		
    	if(mutList==null || mutList.size()<1) return;
    	if(mutList2==null || mutList2.size()<1) return;
    	if(mutList.size()!=mutList2.size()){
    	  System.err.println("Error,Sample num in forward file list is different with the one in reverse num.");
    	  return;
    	}
    	
		List<ArrayList<String>> out=new ArrayList<ArrayList<String>>();
		ArrayList<String> perRow;
		
		//Save average base mut rate as profile
		perRow=new ArrayList<String>();
		perRow.add("Item");	
	    for(SeqMutInfo mut:mutList){	
	      perRow.add(mut.expName.substring(mut.expName.lastIndexOf("/")+1,mut.expName.length()));		
		}
	    out.add(perRow);
		perRow=null;
		
		//========= for forward=========
		//save reads count for each sample
	    perRow=new ArrayList<String>();
	    perRow.add("Forward Reads Count");
	    for(SeqMutInfo mut:mutList){
	      perRow.add(Integer.toString(mut.readsNum));
		}
	    out.add(perRow);
		perRow=null;
		
		//save Average Base Mut
		if(doInterestMut){
		  String mutTag="";
		  for(int m=0;m<interestMut.length;m++){
			  mutTag=mutTag+interestMut[m][0]+">"+interestMut[m][1]+"|";
		  }
		  mutTag=mutTag.substring(0,mutTag.lastIndexOf("|"));
		  perRow=new ArrayList<String>();
		  perRow.add("Forward Average BMF:"+mutTag);
	      for(SeqMutInfo mut:mutList){
	    	 perRow.add(Double.toString(mut.avgTarMutRate));
		  }
	      out.add(perRow);
		  perRow=null;
			
	      perRow=new ArrayList<String>();
	      perRow.add("Forward Average BMF:Other");
	      for(SeqMutInfo mut:mutList){
	    	 perRow.add(Double.toString(mut.avgOtherMutRate));
		  }
	      out.add(perRow);
		  perRow=null;
		}		
		
	    perRow=new ArrayList<String>();
	    perRow.add("Forward Average BMF:All");
	    for(SeqMutInfo mut:mutList){
	      perRow.add(Double.toString(mut.avgMutRate));
		}
	    out.add(perRow);
		perRow=null;
		
		//========= for reverse=========
		//save reads count for each sample
		perRow=new ArrayList<String>();
		perRow.add("Reverse Reads Count");
		for(SeqMutInfo mut2:mutList2){
		   perRow.add(Integer.toString(mut2.readsNum));
		}
		out.add(perRow);
		perRow=null;
				
		//save Average Base Mut
		if(doInterestMut){
		   String mutTag="";
	 	   for(int m=0;m<interestMut.length;m++){
			  mutTag=mutTag+interestMut[m][0]+">"+interestMut[m][1]+"|";
		   }
	 	   mutTag=mutTag.substring(0,mutTag.lastIndexOf("|"));
		   perRow=new ArrayList<String>();
		   perRow.add("Reverse Average BMF:"+mutTag);
		   for(SeqMutInfo mut2:mutList2){
			  perRow.add(Double.toString(mut2.avgTarMutRate));
		   }
		   out.add(perRow);
		   perRow=null;
					
		   perRow=new ArrayList<String>();
		   perRow.add("Reverse Average BMF:Other");
		   for(SeqMutInfo mut2:mutList2){
		   	 perRow.add(Double.toString(mut2.avgOtherMutRate));
		   }
	       out.add(perRow);
		   perRow=null;
		}		
				
		perRow=new ArrayList<String>();
		perRow.add("Reverse Average BMF:All");
		for(SeqMutInfo mut2:mutList2){
		    perRow.add(Double.toString(mut2.avgMutRate));
		}
		out.add(perRow);
		perRow=null;
		
		//========= for forward-reverse=========
		//save Average Base Mut
		double avgMutRate=0.0d;
		if(doInterestMut){
			  String mutTag="";
			  for(int m=0;m<interestMut.length;m++){
				  mutTag=mutTag+interestMut[m][0]+">"+interestMut[m][1]+"|";
			  }
			  mutTag=mutTag.substring(0,mutTag.lastIndexOf("|"));
			  perRow=new ArrayList<String>();
			  perRow.add("Forward-Reverse Average BMF:"+mutTag);
			  avgMutRate=0.0d;
		      for(int i=0;i<mutList.size();i++){
		    	 avgMutRate=mutList.get(i).sumTarMutRate+mutList2.get(i).sumTarMutRate;
		    	 avgMutRate=avgMutRate/(mutList.get(i).tarBaseNum+mutList2.get(i).tarBaseNum);
		    	 perRow.add(Double.toString(avgMutRate));
			  }
		      out.add(perRow);
			  perRow=null;
					
		      perRow=new ArrayList<String>();
		      perRow.add("Forward-Reverse Average BMF:Other");
		      avgMutRate=0.0d;
		      for(int i=0;i<mutList.size();i++){
		    	 avgMutRate=mutList.get(i).sumOtherMutRate+mutList2.get(i).sumOtherMutRate;
		    	 avgMutRate=avgMutRate/(mutList.get(i).otherBaseNum+mutList2.get(i).otherBaseNum);
		    	 perRow.add(Double.toString(avgMutRate));
			  }
		      out.add(perRow);
			  perRow=null;
		}		
				
		perRow=new ArrayList<String>();
		perRow.add("Forward-Reverse Average BMF:All");
		avgMutRate=0.0d;
	    for(int i=0;i<mutList.size();i++){
	    	 avgMutRate=mutList.get(i).sumMutRate+mutList2.get(i).sumMutRate;
	    	 avgMutRate=avgMutRate/(mutList.get(i).baseNum+mutList2.get(i).baseNum);
	    	 perRow.add(Double.toString(avgMutRate));
		}
		out.add(perRow);
		perRow=null;
		
		String outFile="";
		if(outTag==null) outTag="profile";
				
		if(outDir!=null && new File(outDir).exists())
		  outFile=outDir+"/AverageBaseMutFrequency-"+outTag+".txt";
		else
		  outFile="AverageBaseMutFrequency-"+outTag+".txt";
		
		FileOperate.saveMatrixList(out, outFile);
		out=null;
	}
    
    public void saveSeqSNP(List<BaseSNP> snpList, String outDir, String outTag){
		
    	if(snpList==null || snpList.size()<1) return;
		
		List<ArrayList<String>> out=new ArrayList<ArrayList<String>>();
        ArrayList<String> perRow;
        perRow=new ArrayList<String>();
		perRow.add("Base");
		perRow.add("refPos");
		perRow.add("SNP_Rate");
		perRow.add("G");
		perRow.add("C");
		perRow.add("A");
		perRow.add("T");
		perRow.add("Deletion");
		perRow.add("Insertion");
		out.add(perRow);
		for(int i=0;i<snpList.size();i++){
			perRow=new ArrayList<String>();			
			perRow.add(snpList.get(i).refBase);
			perRow.add(Integer.toString(snpList.get(i).refPos));
			perRow.add(Double.toString(snpList.get(i).snpRate));
			perRow.add(Integer.toString(snpList.get(i).altBaseStat.g));
			perRow.add(Integer.toString(snpList.get(i).altBaseStat.c));
			perRow.add(Integer.toString(snpList.get(i).altBaseStat.a));
			perRow.add(Integer.toString(snpList.get(i).altBaseStat.t));
			perRow.add(Integer.toString(snpList.get(i).altBaseStat.del));
			perRow.add(Integer.toString(snpList.get(i).altBaseStat.insert));
			out.add(perRow);
			perRow=null;
		}
		
		String outFile="";
		if(outTag==null) outTag="";
				
		if(outDir!=null)
		  outFile=outDir+"/SeqSNPs-"+outTag+".txt";
		else
		  outFile="SeqSNPs-"+outTag+".txt";
		
		FileOperate.saveMatrixList(out, outFile);
		out=null;

	}
    
}
