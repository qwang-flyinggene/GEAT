package org.geatools.data.load;
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

import org.geatools.data.structure.ChrInfo;
import org.geatools.data.structure.ChrSite;
import org.geatools.operation.FileOperation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class LoadChrInfo {
	
	 static int bandChrCol=0;
	 static int bandStartCol=1;
	 static int bandEndCol=2;
	 static int bandNameCol=3;
	 static int bandTypeCol=4;
	 
	 static int chrInfoNameCol=0;
	 static int chrInfoLengthCol=1;
	 
	 static String chrRegex="[Cc][hH][rR]([1-9]|[1-9][0-9]|[Xx]|[Yy]|[Mm])";
	 static String chrNonXYMRegex="[Cc][hH][rR]([1-9]|[1-9][0-9])";	
	 
	 public static String getChrRegex(){
		 return chrRegex;
	 }
	 
	 public static void setChrRegex(String chromRegx){
		 chrRegex=chromRegx;
	 }
	 
	 public static String getChrNonXYMRegex(){
		 return chrNonXYMRegex;
	 }
	 
	 public static void setChrNonXYMRegex(String chromNonXYMRegx){
		 chrNonXYMRegex=chromNonXYMRegx;
	 }
	 
	 public static List <ChrInfo> getChrInfo(String chrInfoFile){		  
		    
		  List <ChrInfo> chrInfoList=new ArrayList <ChrInfo> ();
		  ChrInfo chrInfo;
		  List <ArrayList <String>> chromInfo0=FileOperation.getMatrixFromFile(chrInfoFile);
		  List<String> chrList=new ArrayList<String>(); 
		  for(int i=0;i<chromInfo0.size(); i++){
			 chrInfo=new ChrInfo ();
			 if(chromInfo0.get(i).get(0).matches(chrRegex)){   
			    chrInfo.name=chromInfo0.get(i).get(chrInfoNameCol).trim();
			    chrInfo.length=Integer.parseInt(chromInfo0.get(i).get(chrInfoLengthCol).trim());
				if(!chrList.contains(chrInfo.name)){
				  chrList.add(chrInfo.name);
				  chrInfoList.add(chrInfo);
				}
			 }	
		   }
		   chrList=null;
		   chrInfoList=setChrNum(chrInfoList);
		  
		   return chrInfoList;
			
	 }
	 
	 public static  List <ChrSite> getChrBands(String cytobandFile){				
		   
		 ArrayList <ChrSite> chrBandsList=new ArrayList<ChrSite>();
		 ChrSite chrBand=null;
		 List<ArrayList<String>> bandList=FileOperation.getMatrixFromFile(cytobandFile);
	         
	     String bandChr="";
	     String bandName="";
	     String bandType="";
	     int bandStart=0;
	     int bandEnd=0;

	     for (int i=0;i<bandList.size();i++) {            
	      	bandChr=bandList.get(i).get(bandChrCol).trim();
	     	if(bandChr.matches(chrRegex)){
		        bandStart=Integer.parseInt(bandList.get(i).get(bandStartCol).trim());
		       	bandEnd=Integer.parseInt(bandList.get(i).get(bandEndCol).trim());
		        bandName=bandList.get(i).get(bandNameCol).trim();
		        bandType=bandList.get(i).get(bandTypeCol).trim(); 
		            
		        chrBand=new ChrSite();
		        chrBand.chr=bandChr;
		        chrBand.name=bandName;
		        chrBand.strand="";
		        chrBand.type=bandType;
		        chrBand.chrStart=bandStart;
		        chrBand.chrEnd=bandEnd;
		        chrBandsList.add(chrBand);	             
	     	}
	     }
	     	
	     return chrBandsList;
			
	 }	

	 public static List<String> getChrNames(List <ChrInfo> chrInfoList){
	     	
	     List<String> chrList=new ArrayList<String>(); 
	     /*
	     String []chrArray=new String[chrInfoList.size()];
	     String chr="";
	     for (int i=0;i<chrInfoList.size();i++) { 
	        chr=chrInfoList.get(i).name;      		 
	        chrArray[chrInfoList.get(i).num-1]=chr;
	     }
	     
	     for(int i=0;i<chrArray.length;i++){
	    	chrList.add(chrArray[i]); 
	     }
	     chrArray=null;
	     */
	     
	     for (int i=0;i<chrInfoList.size();i++) { 
	    	chrList.add(chrInfoList.get(i).name);      		 
		 }
	     
	     return chrList;			
	 }	
	 
	 public static List <ChrInfo> setChrNum(List <ChrInfo> chrInfo){		  
		
		 int chrNonXYMNum=0;
			
		 for(ChrInfo chr:chrInfo){
			if(chr.name.matches(chrNonXYMRegex)) chrNonXYMNum=chrNonXYMNum+1;
		 }
		 
		 String chrName0;
		 String chrName;
		 int chrNum;
		 for(int i=0; i<chrInfo.size();i++){
			chrName=chrInfo.get(i).name;
            chrName0=chrName.substring(3,chrName.length());
            if(chrName0.equalsIgnoreCase("X")) 
            	chrNum=chrNonXYMNum+1;
            else if(chrName0.equalsIgnoreCase("Y")) 
            	chrNum=chrNonXYMNum+2;
            else if(chrName0.equalsIgnoreCase("M")) 
            	chrNum=chrNonXYMNum+3;
            else
            	chrNum=Integer.parseInt(chrName0);
            
            chrInfo.get(i).num=chrNum;
		 }
		 
		 return chrInfo;
	 }
	 
	 public static List <ChrInfo> sortChrByNum(List <ChrInfo> chrInfoList){

		 Collections.sort(chrInfoList, new Comparator<ChrInfo>() {
		     public int compare(ChrInfo  chr1, ChrInfo  chr2){
		        return  Integer.valueOf(chr1.num).compareTo(chr2.num);
		     }
		 });
		 
		 return chrInfoList;
	 }
	 
	 
}
