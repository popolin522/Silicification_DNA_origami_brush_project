# Create input data configuration, with style "full"
import numpy as np
import os
def Write_coordinate(Backbone_length,Chain_length, backbone_molid=1,chain_molid=2,q_rod=-2,q_chain=-2,with_Solvent=False,apply_blockage=True,Config="honeycomb_6HB_atomstyle_full"):
    '''
    Write out the xyz coordinates of 6HB. charged species are added according to the number beads on the 6HB*charges to preserve electroneutrality
    apply_blockage = True, put two bead at the end of 6HB
    '''
    if Config!="honeycomb_6HB_atomstyle_full":
        print("Implementation Error")
    else:
        pass
    i=1 #atomid
    type2=2
    type3=3
    #gf1 = i
    gx = 0
    gy = 0
    gz = 0
    basis = [(1,1), (2,0), (1,-1), (-1,-1), (-2,0), (-1,1)] #y,z 
    atoms_count=0
    backbone_list = []
    chain_list = []
    solvent_list = []
    string_to_write = []
    string_to_write.append("\nAtoms\n\n")

    if apply_blockage==True: #add two beads at the end
        string_to_write.append("%8i %8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,1,5,0,0,0,0)) #type = 5
        i+=1
        string_to_write.append("%8i %8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,1,5,0,-Backbone_length+1,0,0))
        i+=1
        atoms_count+= 2
    else:
        pass

    for face_i in basis:
        gx = 0
        gy = 0
        gz = 0

        x_increment = 1
        y_increment = 1/2
        z_increment = (3**0.5)/2

        y_direction, z_direction = face_i
        y_increment *=y_direction
        z_increment *=z_direction
        for rod_i in range(Backbone_length): #6 faces, honeycomb shape
            gy = y_increment
            gz = z_increment
            string_to_write.append("%8i %8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,backbone_molid,type3,q_rod,gx,gy,gz)   )
            backbone_list.append(i)
            i+=1   
            for graft_i in range(Chain_length):
                gy += y_increment 
                gz += z_increment
                string_to_write.append("%8i %8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,chain_molid,type2,q_chain,gx,gy,gz))
                chain_list.append(i)
                i += 1
            chain_molid += 1
            gx -=x_increment
    if with_Solvent==True:
            solvent_number = 30*30*30
            solvent_size = int(round(solvent_number**(1./3)))
            residue = solvent_number - solvent_size**3
            print(solvent_size)
            print(residue) 
            for x in range(solvent_size+residue):
                for y in range(solvent_size):
                    #for z in range(solvent_size):
                    for  z in range(solvent_size):
                        string_to_write.append("%8i %8i %8i %8i %14.6f %14.6f %14.6f\n" % (i,chain_molid+1,4,1,-40+x,y,-120+z)) #shifting solvent initial condition
                        solvent_list.append(i)
                        i+=1

    atoms_count+=len(backbone_list)+len(chain_list)+len(solvent_list)
    backbone_list = np.array(backbone_list).reshape((6,Backbone_length))
    chain_list = np.array(chain_list).reshape((6,Backbone_length,Chain_length))
    return string_to_write,backbone_list, chain_list,atoms_count

def Write_bond(Backbone_length,Chain_length, backbone_list,chain_list,string_to_write,with_hairpin=True,apply_blockage=True):
    bond_id=1
    bonds_count = 0
    string_to_write.append("\nBonds\n\n")
    #How many bonds to form? 6 X backbone length for inner bond
    for i in range(Backbone_length):
        j=0
        while j < 6:
            if j != 5:
                k= j+1
            else:
                k = 0
            string_to_write.append("{}         {} {} {}\n".format(bond_id, 2,   backbone_list[j][i],  backbone_list[k][i] ))
            bonds_count+=1
            j+=1
            bond_id+=1
        
    for i in range(6):
        j=0
        while j < Backbone_length-1: #number of section-section bonds = Backbone_length -1
            k=j+1

            string_to_write.append("{}         {} {} {}\n".format(bond_id, 2,   backbone_list[i][j],  backbone_list[i][k] ))
            bonds_count+=1
            j+=1
            bond_id+=1 
    #Setting FENE Bonds
    for i in range(6):
        for j in range(Backbone_length):
            k=0
            string_to_write.append("{}         {} {} {}\n".format(bond_id, 1,   backbone_list[i][j],chain_list[i][j][k] )) #bond between backbone and the first bead
            bonds_count+=1
            bond_id+=1
            while k < Chain_length-1: #number of FENE brush bead-brush bead bonds per chain,
                l=k+1    
                string_to_write.append("{}         {} {} {}\n".format(bond_id, 1, chain_list[i][j][k],chain_list[i][j][l] ))
                k+=1
                bond_id+=1
                bonds_count+=1

    if with_hairpin:
        id_that_shouldnt_be_deleted=[]
        for i in range(6):
            for j in range(Backbone_length-1):
                k=0
                while k < Chain_length: #number of FENE brush bead-brush bead bonds per chain,
                    if (k%2) ==0: #if it's an even number bead
                       string_to_write.append("{}         {} {} {}\n".format(bond_id, 1, chain_list[i][j][k],chain_list[i][j+1][k] ))
                       id_that_shouldnt_be_deleted.append(chain_list[i][j][k])
                       id_that_shouldnt_be_deleted.append(chain_list[i][j+1][k])
                       bond_id+=1
                       bonds_count+=1
                    else:
                        pass
                    
                    k+=1
    
    #    
    else:
        id_that_shouldnt_be_deleted = None
        pass

    if apply_blockage==True:
        for i in range(6):
            string_to_write.append("{}         {} {} {}\n".format(bond_id, 1, 1, backbone_list[i][0])) #bonds for the first blockage
            bond_id+=1
            bonds_count+=1
            string_to_write.append("{}         {} {} {}\n".format(bond_id, 1, 2, backbone_list[i][-1])) #bonds for the second blockage
            bond_id+=1
            bonds_count+=1     
    return string_to_write,bonds_count,id_that_shouldnt_be_deleted

def Write_angle(Backbone_length,Chain_length, backbone_list,chain_list,string_to_write,is_dsDNA=False,with_hairpin=True):
    angle_id=1
    string_to_write.append("\nAngles\n\n")
    angles_count = 0
    for i in range(6):
        j=0
        while j < Backbone_length-2: #number of angle bonds = Backbone_length -2
            k=j+1
            l=k+1
            string_to_write.append("{}         {} {} {} {}\n".format(angle_id, 1,   backbone_list[i][j],  backbone_list[i][k],backbone_list[i][l] ))
            j+=1
            angle_id+=1
            angles_count+=1

    if is_dsDNA:
        print("structure: dsDNA")
        #print(chain_list.shape) #it's like [6,backbonelength,chainlength]
        for n in range(6):
            for i in range(Backbone_length):
                j=0
                while j < Chain_length-2: #number of angle bonds per chain = Chain_length -2
                    k=j+1
                    l=k+1
                    string_to_write.append("{}         {} {} {} {}\n".format(angle_id, 1,   chain_list[n][i][j],  chain_list[n][i][k],chain_list[n][i][l] ))
                    j+=1
                    angle_id+=1
                    angles_count+=1
    else:
        print("structure: ssDNA")
    if with_hairpin:
        for n in range(6):
            for i in range(Backbone_length-2):
                j=0
                while j < Chain_length-1: #number of angle bonds per chain = Chain_length -2
                    k=i+1
                    l=k+1
                    string_to_write.append("{}         {} {} {} {}\n".format(angle_id, 1,   chain_list[n][i][j],  chain_list[n][k][j],chain_list[n][l][j] ))
                    j+=1
                    angle_id+=1
                    angles_count+=1
    for i in range(Backbone_length):
        j=0
        while j < 6: #number of angle bonds = Backbone_length -2
            strand_list = ['012','123','234','345','450','501'] #need to make this list circular... this is just a workaround. so 012 means angle potential on strands 0, 1, and 2 
            to_be_applied_angle_potential_idx = strand_list[j]
            k=int(to_be_applied_angle_potential_idx[1])
            l=int(to_be_applied_angle_potential_idx[2])
            string_to_write.append("{}         {} {} {} {}\n".format(angle_id, 2,   backbone_list[j][i],  backbone_list[k][i],backbone_list[l][i] ))
            j+=1
            angle_id+=1
            angles_count+=1

    return  string_to_write, angles_count


def Define_fine_architecture(Backbone_length,id_that_shouldnt_be_deleted=None,is_dumbbell=True,with_hairpin=True):
    '''
    Use this function to define finely coarse-graiened structure
    '''
    delete_molid = ""
    keep_molid_self_and_neighbor=[]
    keep_molid_self_and_neighbor_2 = []
    print("Copy this to lammps input script")
    
    count = 0
    actual_count = 0
    midsection_count = 0
    backbone_iterator = 0
    #for molid in range(2,Backbone_length*6+2):
    molid=2
    if is_dumbbell:

            while molid<Backbone_length*6+2:
                count = backbone_iterator
                midsection_count=0
                actual_count=0
                while actual_count<Backbone_length:
                    if count%6==0:
                        keep_molid_self_and_neighbor.append(molid)
                        pass
                    else:
                        delete_molid += str(molid)
                        delete_molid += " "

                    if midsection_count<30:
                        keep_molid_self_and_neighbor_2.append(molid)
                    elif midsection_count>=Backbone_length-30:
                        keep_molid_self_and_neighbor_2.append(molid)
                        pass
                    else:
                        delete_molid += str(molid)
                        delete_molid += " "

                    count+=1
                    midsection_count+=1
                    molid+=1
                    actual_count+=1

                backbone_iterator+=4
    else:
            while molid<Backbone_length*6+2:
                count = backbone_iterator
                midsection_count=0
                actual_count=0
                while actual_count<Backbone_length:
                    if count%6==0:
                        pass
                    else:
                        delete_molid += str(molid)
                        delete_molid += " "
                    count+=1
                    midsection_count+=1
                    molid+=1
                    actual_count+=1

                backbone_iterator+=4

    if with_hairpin:
        keep_molid = list(set(keep_molid_self_and_neighbor).intersection(keep_molid_self_and_neighbor_2))

        keep_id  = ""
        for id in id_that_shouldnt_be_deleted:
           keep_id+= str(id)
           keep_id += " "
        keep_molid_self_and_neighbor=""
        for mol in keep_molid:
            keep_molid_self_and_neighbor+=str(mol)
            keep_molid_self_and_neighbor += " "
            keep_molid_self_and_neighbor+=str(mol+1)
            keep_molid_self_and_neighbor += " "
            keep_molid_self_and_neighbor+=str(mol+2)
            keep_molid_self_and_neighbor += " "
        
        #Manual adding. The stuff that I added manually
        manual_delete = "195 196 774 775"
        keep_molid_self_and_neighbor+=str(193)
        keep_molid_self_and_neighbor += " "
        keep_molid_self_and_neighbor+=str(192)
        keep_molid_self_and_neighbor += " "
        keep_molid_self_and_neighbor+=str(772)
        keep_molid_self_and_neighbor += " "
        keep_molid_self_and_neighbor+=str(771)
        keep_molid_self_and_neighbor += " "#manually added because it's more easier to do it now

        print("group graft type 2\ngroup keep_mol molecule {}\ngroup delete_chains subtract graft keep_mol\ndelete_atoms group delete_chains bond yes\ngroup keep_id id {}\ngroup delete_mol molecule {}\n group delete_chains subtract delete_mol keep_id\n delete_atoms group delete_chains bond yes\ngroup delete_mol delete\ngroup delete_chains delete\ngroup keep_id delete\n".format(keep_molid_self_and_neighbor,keep_id,delete_molid))
        print("\ngroup manual_delete molecule {}\ndelete_atoms group manual_delete bond yes\ngroup manual_delete delete\n".format(manual_delete))
        #delete side chain in the beginning
        print("group delete_id id 88	604	1120	1636	2152	2668	3184	3700	4216	4732	5248	5764	6280	6796	7312	7828	8344	8860	9376	9892	10408	10924	11440	11956	12472	12988	13504	14020	14536	15052	15568	16084\ngroup delete_id  id 174	690	1206	1722	2238	2754	3270	3786	4302	4818	5334	5850	6366	6882	7398	7914	8430	8946	9462	9978	10494	11010	11526	12042	12558	13074	13590	14106	14622	15138	15654	16170\ngroup delete_id  id 16858	17374	17890	18406	18922	19438	19954	20470	20986	21502	22018	22534	23050	23566	24082	24598	25114	25630	26146	26662	27178	27694	28210	28726	29242	29758	30274	30790	31306	31822	32338	32854		33628	34144	34660	35176	35692	36208	36724	37240	37756	38272	38788	39304	39820	40336	40852	41368	41884	42400	42916	43432	43948	44464	44980	45496	46012	46528	47044	47560	48076	48592	49108	49624\ngroup delete_id  id 16944	17460	17976	18492	19008	19524	20040	20556	21072	21588	22104	22620	23136	23652	24168	24684	25200	25716	26232	26748	27264	27780	28296	28812	29328	29844	30360	30876	31392	31908	32424	32940		33714	34230	34746	35262	35778	36294	36810	37326	37842	38358	38874	39390	39906	40422	40938	41454	41970	42486	43002	43518	44034	44550	45066	45582	46098	46614	47130	47646	48162	48678	49194	49710\ngroup delete_id  id 49882	50398	50914	51430	51946	52462	52978	53494	54010	54526	55042	55558	56074	56590	57106	57622	58138	58654	59170	59686	60202	60718	61234	61750	62266	62782	63298	63814	64330	64846	65362	65878		66652	67168	67684	68200	68716	69232	69748	70264	70780	71296	71812	72328	72844	73360	73876	74392	74908	75424	75940	76456	76972	77488	78004	78520	79036	79552	80068	80584	81100	81616	82132	82648\ngroup delete_id  id 49968	50484	51000	51516	52032	52548	53064	53580	54096	54612	55128	55644	56160	56676	57192	57708	58224	58740	59256	59772	60288	60804	61320	61836	62352	62868	63384	63900	64416	64932	65448	65964		66738	67254	67770	68286	68802	69318	69834	70350	70866	71382	71898	72414	72930	73446	73962	74478	74994	75510	76026	76542	77058	77574	78090	78606	79122	79638	80154	80670	81186	81702	82218	82734\ngroup delete_id  id 83422	83938	84454	84970	85486	86002	86518	87034	87550	88066	88582	89098	89614	90130	90646	91162	91678	92194	92710	93226	93742	94258	94774	95290	95806	96322	96838	97354	97870	98386	98902	99418		16248	16342	66136	66222\ngroup delete_id  id 83508	84024	84540	85056	85572	86088	86604	87120	87636	88152	88668	89184	89700	90216	90732	91248	91764	92280	92796	93312	93828	94344	94860	95376	95892	96408	96924	97440	97956	98472	98988	99504")					

        print('\ndelete_atoms group delete_id bond yes')
        print('\ngroup delete_id delete')
    else:

        print("group delete_mol molecule {}\ndelete_atoms group delete_mol bond yes\ngroup delete_mol delete\n".format(delete_molid))


if __name__ == "__main__":
    Config="honeycomb_6HB_atomstyle_full"
    is_dsDNA = False #change here
    is_dumbbell = True #change here
    with_hairpin = False
    with_Solvent= True
    apply_blockage=True
    if is_dsDNA:
        q_chain=-11.8*5.3 
        file_name = 'dsDNA.data'
    else:
        q_chain = -3.2*5.3 
        file_name = 'ssDNA.data'
    
    if with_hairpin:
        file_name = 'dsDNA_with_hairpin.data'   
    Backbone_length=193 #385nm~192.5beads
    Chain_length=85
    string_to_write, backbone_list,chain_list,atoms_count = Write_coordinate(Backbone_length,Chain_length,with_Solvent=with_Solvent,q_chain=q_chain,q_rod=-11.8*5.3,apply_blockage=apply_blockage) #unit testing , q_chain=-20 for dsdna, and q_chain=-10 for ssdna
    string_to_write,bonds_count,id_that_shouldnt_be_deleted = Write_bond(Backbone_length,Chain_length,backbone_list,chain_list,string_to_write,with_hairpin=with_hairpin,apply_blockage=apply_blockage)
    string_to_write,angles_count = Write_angle(Backbone_length,Chain_length, backbone_list,chain_list,string_to_write,is_dsDNA=is_dsDNA,with_hairpin=with_hairpin)
    Define_fine_architecture(Backbone_length,id_that_shouldnt_be_deleted=id_that_shouldnt_be_deleted, is_dumbbell=is_dumbbell,with_hairpin=with_hairpin)

    with open(file_name, 'w') as f:
        f.write("LAMMPS data file for {}\n".format(Config))
        f.write("{} atoms\n".format(atoms_count))
        f.write("6 atom types\n")
        f.write("{} bonds\n".format(bonds_count))
        f.write("2 bond types\n")
        f.write("{} angles\n".format(angles_count))
        #f.write("1 angle types\n")
        f.write("2 angle types\n") #add addditional angle type 2 to maintain the honeycomb structure
        f.write("-250 150 xlo xhi\n-200 200 ylo yhi\n-200 200 zlo zhi\n")
        f.write("\nMasses\n\n1 1\n2 1\n3 1\n4 1\n5 1\n6 1\n")
        for line in string_to_write:
            f.write(f"{line}")