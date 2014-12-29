#!/bin/sh

############################################################################
#
# Template script for mirroring PDB SNAPSHOTS FTP archive using rsync
#
############################################################################

# This script is being provided to PDB users as a template for using rsync
# to mirror the Snapshots FTP archive from an anonymous rsync server. You may want
# to review rsync documentation for options that better suit your needs.
#
#
#                     *******************************
#
#
# Please read on below and select the rsync option that best suits your
# requirements.  Depending on your network connection rsync'ing the entire
# archive can take a long time (24+ hours). You may only want to make a
# copy of a portion of the archive and some of the most popular sections
# of the archive are highlighted below.
#
#
#                     *********************************
# Last modified:	December 3, 2014

# Comment out the following line if you no longer want to see this info
echo "Prior to first use, you must edit this script to choose a server name, port number, and rsync option!"

# Change the following three lines to reflect your local settings.

MIRRORDIR=./               # your top level rsync directory
LOGFILE=rsync-ftp.log      # file for storing logs
RSYNC=/usr/bin/rsync       # location of local rsync

##########################################################################################
#
#        YOU MUST UNCOMMENT YOUR CHOICE OF SERVER AND CORRESPONDING PORT BELOW
#
SERVER=snapshotrsync.rcsb.org                       # RCSB server name
PORT=8730                                           # port RCSB server is using
#
#SERVER=pdbjsnap.protein.osaka-u.ac.jp               # PDBj server name
#PORT=873                                            # port PDBj server is using
#
##########################################################################################
#
#        YOU MUST UNCOMMENT THE RYSNC OPTION BELOW THAT MEETS YOUR NEEDS!
#
# Rsync the PDB format coordinate data from FTP snapshot of the PDB on Dec 3, 2014
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT \
#$SERVER::20141203/pub/pdb/data/structures/divided/pdb/ \
#$MIRRORDIR > $LOGFILE 2>/dev/null

#
# Rsync the XML format coordinate data from FTP snapshot of the PDB on Dec 3, 2014
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT \
#$SERVER::20141203/pub/pdb/data/structures/divided/XML/ \
#$MIRRORDIR > $LOGFILE 2>/dev/null

#
# Rsync the mmCIF format coordinate data from FTP snapshot of the PDB on Dec 3, 2014
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT \
#$SERVER::20141203/pub/pdb/data/structures/divided/mmCIF/ \
#$MIRRORDIR > $LOGFILE 2>/dev/null


#                TAKE NOTE  - THESE LATTER OPTIONS WILL TAKE A WHILE --


# Rsync the entire FTP snapshots site (All Snapshots - Approx 1.8TB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::all/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Dec 3, 2014 (Approx 434GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20141203/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 2, 2014 (Approx 292GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20140102/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 1, 2013 (Approx 219GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20130101/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 2, 2012 (Approx 135GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20120102/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jul 7, 2011 (Approx 127GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20110707/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 3, 2011 (Approx 117GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20110103/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 4, 2010 (Approx 98GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20100104/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Mar 16, 2009 (Approx 84GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20090316/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 5, 2009 (Approx 81GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20090105/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 7, 2008 (Approx 71GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20080107/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jul 31, 2007 (Approx 66GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20070731/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 2, 2007 (Approx 59GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20070102/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the entire FTP snapshot of the PDB on Jan 3, 2006 (Approx 46GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20060103/ $MIRRORDIR > $LOGFILE 2>/dev/null

# Rsync the FTP snapshot of the PDB on Jan 6, 2005    (Approx 24GB)
#
#${RSYNC} -rlpt -v -z --delete --port=$PORT $SERVER::20050106/ $MIRRORDIR > $LOGFILE 2>/dev/null

