#!/usr/bin/perl
use strict;
use warnings;

#quality hash

my %quality = ("!",0,
             "\"",1,
             "\#",2,
             "\$",3,
             "%",4,
             "&",5,
             "\'",6,
             "(",7,
             ")",8,
             "*",9,
             "+",10,
             ",",11,
             "-",12,
              ".",13,
              "/",14,
              "0",15,
              "1",16,
              "2",17,
              "3",18,
              "4",19,
              "5",20,
              "6",21,
              "7",22,
              "8",23,
              "9",24,
              ":",25,
              ";",26,
              "<",27,
              "=",28,
              ">",29,
              "?",30,
              "@",31,
              "A",32,
              "B",33,
              "C",34,
              "D",35,
              "E",36,
              "F",37,
              "G",38,
              "H",39,
              "I",40,
              "J",41,
              "K",42,
              );

my $i;
my $count=0;
while (<>) {

       $count++;
       #print "line", $count, "\n";
      chomp;
      my @passed =();
      my $tmp = "";
      my @fields=split(/\t/,$_);
      my $depth = $fields[3];
      my $ref = $fields[2];
      my $bases = $fields[4];
      my $quals = $fields[5];
      $bases=~s/[\+\-][1-9]{1,}[ACTGNactgn]{1}//g;
      $bases=~s/\^.//g;
      #$bases=~s/\*//g;
      $bases=~s/[\[\]><\$]//g;
      $bases=~s/[\+\-][1-9]{1,}[ACTGNactgn]{1}//g;
      for($i=0;$i<length($quals);$i++)
      {
          if($quality{substr($quals,$i,1)})
          {

               if(($quality{substr($quals,$i,1)}>=30) && (substr($bases,$i,1) ne "*"))  ####Checar que la calidad 
               {
                       
                       push(@passed,substr($bases,$i,1));
               }
          }else{ print next;} #print "bad qualities $quals\n";
       }
       if(scalar(@passed) == 0)
       {
          next;
       }
       if(scalar(@passed)==1)
       {
               if($passed[0] eq "." || $passed[0] eq ",")
               {
                       print "$_\t $ref\n";
               }
               else
               {
                       print "$_\t $passed[0]\n";
               }
       }
       else
       {
               my $random_allele =$passed[int(rand(scalar(@passed)-1))];

       #       print "random allele = $random_allele \n";

               if($random_allele eq "." || $random_allele eq ",")
               {
                       print "$_\t $ref\n";
               }

               else
               {
                        print "$_\t $random_allele\n";
               }
       }

}


exit;

