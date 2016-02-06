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

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;
import java.util.ArrayList; 

import org.geatools.data.structure.BlastInfo;
import org.geatools.data.structure.ChrSite;
import org.geatools.data.structure.SeqChrAlignSite;
import org.geatools.data.structure.SeqInfo;
import org.geatools.operation.FileOperate;
import org.geatools.seqprocess.SeqOperation;

public class SeqOffTarget{

	 //int siteDownStreamSize=6;
	 //int siteUpStreamSize=0;
	    
	 static String tmpDir="tmp";
	 
	 public SeqOffTarget(){
	 
	 }
	 
	 public void setTmpDir(String dir){
		tmpDir=dir;
	 }
	 public void delTmpDir(String dir){
		FileOperate.delFolder(dir);
	 }	 
	 
	 public List<ChrSite> getOffTargetSite(String targetSeqFile,String refGenomeFile, 
			 float minIdentity,int minAlignLen,int maxMismatchNum,int maxGapNum){
	   
	    List<ChrSite> offTargetSites=new ArrayList<ChrSite> ();
		ChrSite perSite;
		String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String blastOutFile=tmpDir+"/"+refGenomeFile.substring(
		         refGenomeFile.lastIndexOf("/")+1,refGenomeFile.lastIndexOf("."))
		         +"."+timeStamp+".blastOut";
		
		String blastCMD="blastn -task blastn-short -word_size=7 -evalue 10000";	
	    blastCMD=blastCMD+" -query "+targetSeqFile+" -subject "+refGenomeFile
	    		+" -out "+blastOutFile+" -outfmt 6";	
	    SeqOperation.runBLAST2Seq(blastCMD);
		
		String offTargetChr="";	
		String offTargetStrand="";	
		String offTargetName="OffTarget";
		
		int maxPos=0;
		int minPos=0;
		int tmp1=0;
		int tmp2=0;
		float identity=80.00f;
		int alignLen=12;
		int mismatchNum=2;
		int gapNum=2;
		List<SeqInfo> targetSeqInfo=SeqOperation.getFASTASeqObj(targetSeqFile);
		int targetSeqLen=targetSeqInfo.get(0).seq.length();
		targetSeqInfo=null;
		
		BlastInfo blastInfo=new BlastInfo(); 
	    List <ArrayList <String>> blastOut=FileOperate.getMatrixFromFile(blastOutFile);
		for(int n=0;n<blastOut.size();n++){
	     
		 offTargetChr=blastOut.get(n).get(blastInfo.colSName);
		 identity=Float.parseFloat(blastOut.get(n).get(blastInfo.colIdentity));
		 alignLen=Integer.parseInt(blastOut.get(n).get(blastInfo.colAlignLen));
		 mismatchNum=Integer.parseInt(blastOut.get(n).get(blastInfo.colMismatchNum));
		 gapNum=Integer.parseInt(blastOut.get(n).get(blastInfo.colGapNum));		
		 
		 offTargetStrand="+";
	     tmp1=Integer.parseInt(blastOut.get(n).get(blastInfo.colSStart));
		 tmp2=Integer.parseInt(blastOut.get(n).get(blastInfo.colSEnd));	
	     minPos=tmp1;
		 maxPos=tmp2;
		 if(tmp2<tmp1) {
		   minPos=tmp2;
		   maxPos=tmp1;
		   offTargetStrand="-";
		 }	 
		 if(identity>=minIdentity && alignLen>=minAlignLen && mismatchNum<=maxMismatchNum
				 && gapNum<=maxGapNum){
		    perSite=new ChrSite ();
		    
		    offTargetName="OffTarget";
		    if(identity==100.0f && alignLen==targetSeqLen && mismatchNum==0 && gapNum==0)
		    	offTargetName="Target";
		    
		    perSite.chr=offTargetChr;
		    perSite.name=offTargetName;
		    perSite.strand=offTargetStrand;
		    perSite.chrStart=minPos;
		    perSite.chrEnd=maxPos;
		    offTargetSites.add(perSite);
		    perSite=null;
		  
		 }
		}
		
		blastOut=null;
		
		return offTargetSites;
	 
	  }
	  
	  public List<SeqChrAlignSite> getOffTargetSite(String targetSeqFile,String refGenomeFile, 
			 float minIdentity,float minAlignLenRatio,float maxMismatchRatio,
			 int maxGapNum){
	   
	    List<SeqChrAlignSite> offTargetSite=new ArrayList<SeqChrAlignSite> ();
	    SeqChrAlignSite perSite;
	    String timeStamp=new SimpleDateFormat("yyyyMMdd_HHmmss").format(Calendar.getInstance().getTime());
		String blastOutFile=tmpDir+"/"+refGenomeFile.substring(
		    refGenomeFile.lastIndexOf("/"),refGenomeFile.lastIndexOf("."))
		    +"."+timeStamp+".blastOut";
		String blastCMD="blastn -task blastn-short -word_size=7 -evalue 10000";	
	    blastCMD=blastCMD+" -query "+targetSeqFile+" -subject "+refGenomeFile
	    		+" -out "+blastOutFile+" -outfmt 6";	
	    SeqOperation.runBLAST2Seq(blastCMD);
		
		String offTargetChr="";	
		String offTargetStrand="";
		String offTargetName="OffTarget";
		
		List<SeqInfo> targetSeqInfo=SeqOperation.getFASTASeqObj(targetSeqFile);
		int targetSeqLen=targetSeqInfo.get(0).seq.length();
		targetSeqInfo=null;
		int minAlignLen=(int) (minAlignLenRatio*targetSeqLen);
		int maxMismatchNum=2;	
		
		int maxPos=0;
		int minPos=0;
		int tmp1=0;
		int tmp2=0;
		float identity=80.0f;
		int alignLen=12;
		int mismatchNum=2;
		int gapNum=2;
		float bitScore=0.0f;
		BlastInfo blastInfo=new BlastInfo(); 
	    List <ArrayList <String>> blastOut=FileOperate.getMatrixFromFile(blastOutFile);
		for(int n=0;n<blastOut.size();n++){
	     
			 offTargetChr=blastOut.get(n).get(blastInfo.colSName);
			 identity=Float.parseFloat(blastOut.get(n).get(blastInfo.colIdentity));
			 alignLen=Integer.parseInt(blastOut.get(n).get(blastInfo.colAlignLen));
			 mismatchNum=Integer.parseInt(blastOut.get(n).get(blastInfo.colMismatchNum));
			 maxMismatchNum=(int) Math.ceil(alignLen*maxMismatchRatio);
			 gapNum=Integer.parseInt(blastOut.get(n).get(blastInfo.colGapNum));
			 bitScore=Float.parseFloat(blastOut.get(n).get(blastInfo.colBitScore));
					 
			 offTargetStrand="+";
		     tmp1=Integer.parseInt(blastOut.get(n).get(blastInfo.colSStart));
			 tmp2=Integer.parseInt(blastOut.get(n).get(blastInfo.colSEnd));	
		     minPos=tmp1;
			 maxPos=tmp2;
			 if(tmp2<tmp1){
			   minPos=tmp2;
			   maxPos=tmp1;
			   offTargetStrand="-";
			 }	 
			 if(identity>=minIdentity && alignLen>=minAlignLen && mismatchNum<=maxMismatchNum
					 && gapNum<=maxGapNum){
			    
				perSite=new SeqChrAlignSite ();
			    
			    offTargetName="OffTarget";
			    if(identity==100.0f && alignLen==targetSeqLen && mismatchNum==0 && gapNum==0)
			    	offTargetName="Target";
			    
			    perSite.chr=offTargetChr;
			    perSite.name=offTargetName;
			    perSite.strand=offTargetStrand;
			    perSite.chrStart=minPos;
			    perSite.chrEnd=maxPos;
			    perSite.mismatchNum=mismatchNum;
			    perSite.alignLen=alignLen;
			    perSite.score=bitScore;
			    perSite.identity=identity;
			    
			    offTargetSite.add(perSite);
			    perSite=null;
			  
			 }
		}
		
		blastOut=null;
		
		return offTargetSite;
	 
	  }
	  
	  public void sortSiteByScore(List <SeqChrAlignSite> siteList,boolean desc){
	
		 Collections.sort(siteList, new SeqChrAlignSite.CompScore(desc)); 
	
	  }
	  
	  public void sortSiteByMismatchNum(List <SeqChrAlignSite> siteList,boolean desc){
	
		 Collections.sort(siteList, new SeqChrAlignSite.CompMismatch(desc)); 
	
	  }
	  
	 
 }