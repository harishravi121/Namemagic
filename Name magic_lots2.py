import random
import time
start=time.time()
w='ATCG'
flag=0
NAME="DR HARISH MCCALTEO RAVI"

def translate(seq): 
       
    table = { 
        'A':'a=F/m','C':'c^4m^2+p^2c^2=E^2','H':'h=E/f', 'M':'MmG/r^2=F','T':'T=PV/nR','R':'r=sqrt(F/GMm)','D':'D=C/pi','W':'W=mg','I':'i=ln(-1)/pi',
        'S':'S=ut+1/2at^2','V':'V^2-u^2=2aS','U':'u=s/t-.5at',' ':'  ','E':'E=mc^2','L':'L=L0sqrt(1-v^2/c^2)','K':'k=(E-Ef)/T/ln(1/n-1)','N':'n=1/(e^(E-Ef)/kT-1)',
        'B':'B=sqrt((A+B)^2-A^2-2AB)','J':'J=I/A','G':'G=Fr^2/M/m','Y':'Y=cuberoot((X+Y)^3-X^3-3XY^2-3X^2Y)','Z':'R+j(wl-1/wc)','O':'O(nlog(n))','F':'F=96500 C/mol',
        'P':'P=D-e0E','.':'.','X':'X=cuberoot((X+Y)^3-Y^3-3XY^2-3X^2Y)'
    } 
    physics ="" 
   
    for i in range(0, len(seq)): 
        
        physics+= table[seq[i]]+' ' 
    return physics 


physics=NAME
p2=translate(physics)
print(physics)
print(p2)
print()
import random
import time
start=time.time()
w='ATCG'
flag=0
def translate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            codon = seq[i:i + 3] 
            protein+= table[codon] 
    return protein 

def revtranslate(seq): 
       
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        'XXX':' ', 'ATX':'U','AXT':'O','XAG':'X','XTG':'B',
        'XTA':'J','XTA':'J','XGA':'Z','.':'.'
    } 
    DNA ="" 
    table2 = {v: k for k, v in table.items()}
    for i in range(0, len(seq)): 
         
        DNA+= table2[seq[i]] 
    return DNA
def RNA(seq): 
       
    table3 = { 
        'A':'U','T':'A','C':'G','G':'C','X':'X','.':'.'
    } 
    RNA ="" 
    
    for i in range(0, len(seq)): 
         
        RNA+= table3[seq[i]] 
    return RNA


protein=NAME
print('Protein ',protein)
DNA=revtranslate(protein)
RNA=RNA(DNA)
print('DNA ',DNA)
print('RNA ',RNA)

import random
import time
start=time.time()
flag=0
#Array of words


def convert(lst): 
    return (lst[0].split()) 
  
# Driver code 
lst =  [NAME] 
print()
a= convert(lst)
def sentancify(s): #Given an array of characters, it makes the word
    s2=''
    for i in range(len(s)):
        s2=s2+s[i]+' '
    return s2

def wordify(w): #Given an array of characters, it makes the word
    c=''
    for i in range(len(w)):
        c=c+w[i]
    return c


for i in range(2):
    b=random.sample(a,len(a)) #sample the word array randomly
    #print(sentancify(b))
    for i in range(2):
        s3=''
        for i in range(len(a)): #5 Words scrambled
       
            word=b[i]
            c=random.sample(word,len(word)) #sample the word array randomly
            f=wordify(c) #It makes the word
            s3=s3+f
            s3=s3+' '

        print(s3)
        


def translate(seq): 
       
    table = { 
        'E':'€','P':'£','Y':'¥','S':'$','N':'₦','R':'₹','B':'₿','M':'ɱ','W':'₩','T':'৳',
        'G':'₾','A':'A','C':'¢','D':'D','F':'F','H':'H','I':'I','J':'J','K':'K','N':'N','O':'O',
        'Q':'Q','U':'U','V':'V','X':'X','Z':'Z',' ':' ','L':'L','.':'.',',':','
    } 
    physics ="" 
   
    for i in range(0, len(seq)): 
        
        physics+= table[seq[i]] 
    return physics 


curr=NAME
p2=translate(curr)

print(p2)
print()
#print(curr)


lst =  [p2] 
print()
a= convert(lst)
def sentancify(s): #Given an array of characters, it makes the word
    s2=''
    for i in range(len(s)):
        s2=s2+s[i]+' '
    return s2

def wordify(w): #Given an array of characters, it makes the word
    c=''
    for i in range(len(w)):
        c=c+w[i]
    return c


for i in range(2):
    b=random.sample(a,len(a)) #sample the word array randomly
    #print(sentancify(b))
    for i in range(2):
        s3=''
        for i in range(len(a)): #5 Words scrambled
       
            word=b[i]
            c=random.sample(word,len(word)) #sample the word array randomly
            f=wordify(c) #It makes the word
            s3=s3+f
            s3=s3+' '

        print(s3)
        
# Python program to implement Morse Code Translator
  
'''
VARIABLE KEY
'cipher' -> 'stores the morse translated form of the english string'
'decipher' -> 'stores the english translated form of the morse string'
'citext' -> 'stores morse code of a single character'
'i' -> 'keeps count of the spaces between morse characters'
'message' -> 'stores the string to be encoded or decoded'
'''
  
# Dictionary representing the morse code chart
MORSE_CODE_DICT = { 'A':'.-', 'B':'-...',
                    'C':'-.-.', 'D':'-..', 'E':'.',
                    'F':'..-.', 'G':'--.', 'H':'....',
                    'I':'..', 'J':'.---', 'K':'-.-',
                    'L':'.-..', 'M':'--', 'N':'-.',
                    'O':'---', 'P':'.--.', 'Q':'--.-',
                    'R':'.-.', 'S':'...', 'T':'-',
                    'U':'..-', 'V':'...-', 'W':'.--',
                    'X':'-..-', 'Y':'-.--', 'Z':'--..',
                    '1':'.----', '2':'..---', '3':'...--',
                    '4':'....-', '5':'.....', '6':'-....',
                    '7':'--...', '8':'---..', '9':'----.',
                    '0':'-----', ', ':'--..--', '.':'.-.-.-',
                    '?':'..--..', '/':'-..-.', '-':'-....-',
                    '(':'-.--.', ')':'-.--.-'}
  
# Function to encrypt the string
# according to the morse code chart
def encrypt(message):
    cipher = ''
    for letter in message:
        if letter != ' ':
  
            # Looks up the dictionary and adds the
            # correspponding morse code
            # along with a space to separate
            # morse codes for different characters
            cipher += MORSE_CODE_DICT[letter] + ' '
        else:
            # 1 space indicates different characters
            # and 2 indicates different words
            cipher += ' '
  
    return cipher
  
# Function to decrypt the string
# from morse to english
def decrypt(message):
  
    # extra space added at the end to access the
    # last morse code
    message += ' '
  
    decipher = ''
    citext = ''
    for letter in message:
  
        # checks for space
        if (letter != ' '):
  
            # counter to keep track of space
            i = 0
  
            # storing morse code of a single character
            citext += letter
  
        # in case of space
        else:
            # if i = 1 that indicates a new character
            i += 1
  
            # if i = 2 that indicates a new word
            if i == 2 :
  
                 # adding space to separate words
                decipher += ' '
            else:
  
                # accessing the keys using their values (reverse of encryption)
                decipher += list(MORSE_CODE_DICT.keys())[list(MORSE_CODE_DICT
                .values()).index(citext)]
                citext = ''
  
    return decipher
  
# Hard-coded driver function to run the program

    
message = NAME
#print(message)
print()
result = encrypt(message.upper())
print (result)

import numpy as np
nums=[0 for i in range(len(NAME))]
nums2=[0 for i in range(len(NAME))]
for i in range(len(nums)):
    
    nums[i]= ord(NAME[i])-64
    nums2[i]=( ord(NAME[i])-64)**2
    
print(nums)
print(nums2) 
