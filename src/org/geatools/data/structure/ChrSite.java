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
package org.geatools.data.structure;

import java.util.Comparator;

public class ChrSite {
	
	public int index;
	public String chr;
	public String name;
	public String type="-";
	public String strand;
	public int chrStart;
	public int chrEnd;
	public String associatedGen="-";
	public String genRelation="-";
	public Double score=null;
	public String seq="";
	
	public static class CompSiteCenter implements Comparator<ChrSite> {
		
		private int mod = 1; //default ascend
		
		public CompSiteCenter() {
			mod = 1; //ascend
		}

		public CompSiteCenter(boolean desc) {
		  if (desc) mod =-1; //descend
		}
		        
		public int compare(ChrSite  site1, ChrSite  site2){
		  return  mod*Integer.valueOf((site1.chrStart+site1.chrEnd)/2)
				  .compareTo((site2.chrStart+site2.chrEnd)/2);
		}
		
	}
	
	public static class CompSiteStart implements Comparator<ChrSite> {
		
		private int mod = 1; //default ascend
		
		public CompSiteStart() {
			mod = 1; //ascend
		}

		public CompSiteStart(boolean desc) {
		  if (desc) mod =-1; //descend
		}
		        
		public int compare(ChrSite  site1, ChrSite  site2){
		  return  mod*Integer.valueOf(site1.chrStart).compareTo(site2.chrStart);
		}
		
	}
	
	public static class CompSiteEnd implements Comparator<ChrSite> {
		
		private int mod = 1; //default ascend
		
		public CompSiteEnd() {
			mod = 1; //ascend
		}

		public CompSiteEnd(boolean desc) {
		  if (desc) mod =-1; //descend
		}
		        
		public int compare(ChrSite  site1, ChrSite  site2){
		  return  mod*Integer.valueOf(site1.chrEnd).compareTo(site2.chrEnd);
		}
		
	}

	
}
