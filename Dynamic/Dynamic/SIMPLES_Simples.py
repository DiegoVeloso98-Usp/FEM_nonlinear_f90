import numpy as np

fnodes = np.array([# cx, cy
                 #[ 2.0, 5.0],
                 #[ 8.0, 5.0]
                 [ 0.0, 0.0],
                 [1.0, 0.0]
                #  [ 0.0,10.0],
                #  [10.0,10.0]
                ])

fconec = np.array([#nó 1 nó 2 properties
                  [0,1,0]
                #   [0,2,0],
                #   [0,3,0],
                #   [2,3,0],
                #   [1,3,0]
                 ])

fproperties = np.array([#young, area
                       [00000.0, 0.0012]
                      ])

nodes = np.array([# cx, cy
                 [ 0.0, 0.0],
                 [1.0, 0.0],
                 [ 1.0,1.0],
                [0.0,1.0]
                ])

conec = np.array([#nó 1..... nó n,properties
                  [0,1,2,0],
                 [0,2,3,0]
                 ])

properties = np.array([#young, ni, h
                       [1000.0, 0.3, 0.1]
                      ])

vinc = np.array([#nó dir val
                 [0,0,0.0],
                 [0,1,0.0],
                 [1,0,0.0],
                 [1,1,0.0],
                ])

nload = np.array([#
                #   [2,1,10.0],
                  [2,1,-0.1]
                 ])

eload = np.array([#elem bx10 bx20 bx11 bx21 bx12 bx22
                  [0,0,0,0,0,0,0],
                  [0,0,0,0,0,0,0]
                 ])

def fforma(ksi,eta):
    # fi = np.array([ksi,eta,1 - ksi - eta])
    # dfidksi = np.array([1,0,-1])
    # dfideta = np.array([0,1,-1])
    fi = np.array([1 - ksi - eta,ksi,eta])
    dfidksi = np.array([-1,1,0])
    dfideta = np.array([-1,0,1])
    return fi,dfidksi,dfideta

#número de passos de carga
npc = 1
#tolerância
tol = 0.000001
#número de nós das fibras
nnodesf = len(fnodes)
#número de fibras
nfibras = len(fconec)
#número de nós da estrutura
nnodes = len(nodes)
#número de elementos da estrutura
nelems = len(conec)
#número de nós por elemento
npe = 3
#número de pontos de integração de Hammer
nph = 12
#número de graus de liberdade
ngl = 2

#identificação das posições das fibras
ifconec = []
for ino in fnodes:
    x1p = ino[0]
    x2p = ino[1]
    for iel in range(nelems):
        x11 = nodes[conec[iel,0],0]
        x21 = nodes[conec[iel,0],1]
        x12 = nodes[conec[iel,1],0]
        x22 = nodes[conec[iel,1],1]
        x13 = nodes[conec[iel,2],0]
        x23 = nodes[conec[iel,2],1]
        k1 = -((-(x13*x22) + x1p*x22 + x12*x23 - x1p*x23 - x12*x2p + x13*x2p)/(x12*x21 - x13*x21 - x11*x22 + x13*x22 + x11*x23 - x12*x23))
        k2 = -((x13*x21 - x1p*x21 - x11*x23 + x1p*x23 + x11*x2p - x13*x2p)/(x12*x21 - x13*x21 - x11*x22 + x13*x22 + x11*x23 - x12*x23))
        k3 = 1 - k1 - k2
        if (k1 >= 0.0) and (k1 <= 1.0) and (k2 >= 0.0) and (k2 <= 1.0) and (k3 >= 0.0) and (k3 <= 1.0):
            ifconec.append([iel,k1,k2])
            break
# print ("ifconec vale: ", ifconec)        


mbx = np.zeros((nelems,ngl*npe))
for iload in eload:
    elem = int(iload[0])
    vbx = iload[1:]
    mbx[elem] = vbx

hammer = {
            "ksi":[0.501426509658179, 0.249286745170910, 0.249286745170910,
                    0.873821971016996, 0.063089014491502, 0.063089014491502,
                    0.053145049844816, 0.310352451033785, 0.636502499121399,
                    0.310352451033785, 0.636502499121399, 0.053145049844816],
            "eta":[0.249286745170910, 0.249286745170910, 0.501426509658179,
                    0.063089014491502, 0.063089014491502, 0.873821971016996,
                    0.310352451033785, 0.636502499121399, 0.053145049844816,
                    0.053145049844816, 0.310352451033785, 0.636502499121399],
            "peso":[0.116786275726379/2.0, 0.116786275726379/2.0, 0.116786275726379/2.0,
                     0.050844906370207/2.0, 0.050844906370207/2.0, 0.050844906370207/2.0,
                     0.082851075618374/2.0, 0.082851075618374/2.0, 0.082851075618374/2.0,
                     0.082851075618374/2.0, 0.082851075618374/2.0, 0.082851075618374/2.0]
         }

y = np.array(nodes)
yf = np.array(fnodes)

iconec = np.zeros(nnodes)
for i in conec:
    for j in range(npe):
        iconec[i[j]] += 1

# print ("iconec vale: ", iconec)
'''
ksieta = np.array([
                   [1.0,0.0],
                   [0.0,1.0],
                   [0.0,0.0]
                  ])
'''
ksieta = np.array([
                   [0.0,0.0],
                   [1.0,0.0],
                   [0.0,1.0]
                  ])
saida = []
saida.append('arquivo no formato acadview \n')
saida.append('nnodes nelems nlistas \n')
saida.append('#\n')
saida.append(f'{nnodes} {nelems} {int(5*npc)}\n')
saida.append('cx cy cz dx dy dz \n')
saida.append('#\n')
for i in range(nnodes):
    saida.append(f'{nodes[i][0]} {nodes[i][1]} 0 0 0 0\n')
saida.append('tipo gaprox no_1....no_npe group\n')
saida.append('#\n')
for iel in range(nelems):
    saida.append(f'2 1 {conec[iel][0]+1} {conec[iel][1]+1} {conec[iel][2]+1} 0\n')
saida.append('listas de resultados\n')
saida.append('nome da lista\n')
saida.append('dx dy dz valor_cores\n')
counter=0
for ipc in range(npc):
    print("passo: ",ipc)
    #variação do resultado
    erro = 1.0
    fext = np.zeros(ngl*nnodes)
    for iload in nload:
        node = int(iload[0])
        dir = int(iload[1])
        val = iload[2]
        fext[ngl * node + dir] += (ipc + 1) / npc * val
        #print ("fext vale: ", fext)
        
    
    for idispla in vinc:
        node = int(idispla[0])
        dir = int(idispla[1])
        val = idispla[2]
        y[node,dir] = nodes[node,dir] + val * (ipc + 1) / npc
    #print ("y vale: ", y)
    
    iter = 0
    while erro > tol:
        fint = np.zeros(ngl*nnodes)
        hessiana = np.zeros((ngl*nnodes,ngl*nnodes))
        for iel in range(2):
            young = properties[conec[iel,-1],0]
            ni = properties[conec[iel,-1],1]
            h = properties[conec[iel,-1],2]

            for ih in range(nph):
                efext = np.zeros(ngl*npe)
                ef = np.zeros(ngl*npe)
                eh = np.zeros((ngl*npe,ngl*npe))

                ksi = hammer["ksi"][ih]
                eta = hammer["eta"][ih]
                peso = hammer["peso"][ih]

                # print(f"Ponto: {ih} ksi: {ksi} eta: {eta} peso: {peso}")

                fi,dfidksi,dfideta = fforma(ksi,eta)
                #print("fi vale:", fi)
                mfi = np.dot(fi.reshape(-1,1),fi.reshape(1,-1))
                #print("mfi vale:", mfi)
                mfifi = np.zeros((ngl*npe,ngl*npe))
                
               
                for ii in range(npe):
                    for jj in range(npe):
                        for kk in range(ngl):
                            mfifi[ngl*ii + kk,ngl*jj + kk] = mfi[ii,jj]

                efext = np.dot(mfifi,(ipc + 1) / npc*mbx[iel])

                A0 = np.zeros((2,2)) 
                A1 = np.zeros((2,2)) 
                for ino in range(npe):
                    A0[0,0] += dfidksi[ino] * nodes[conec[iel,ino],0]
                    A0[0,1] += dfideta[ino] * nodes[conec[iel,ino],0]
                    A0[1,0] += dfidksi[ino] * nodes[conec[iel,ino],1]
                    A0[1,1] += dfideta[ino] * nodes[conec[iel,ino],1]
                    A1[0,0] += dfidksi[ino] * y[conec[iel,ino],0]
                    A1[0,1] += dfideta[ino] * y[conec[iel,ino],0]
                    A1[1,0] += dfidksi[ino] * y[conec[iel,ino],1]
                    A1[1,1] += dfideta[ino] * y[conec[iel,ino],1]


                # print("A0 vale: ", A0)
                # print("A1 vale: ", A1)
                jota = np.linalg.det(A0)
                # print("jota vale: ", jota)
                # exit()
                # print("A0 vale: ", A0)
                # exit()
                #EPT
                A = np.dot(A1,np.linalg.inv(A0))
                C = np.dot(np.transpose(A),A)
                E = 1/2 * (C - np.identity(2))
                # print("E vale: ", E)
                # exit()
                S = young / (1 - ni**2) * np.array([
                                                    [E[0,0] + ni * E[1,1],    (1 - ni) * E[0,1]],
                                                    [   (1 - ni) * E[1,0], ni * E[0,0] + E[1,1]]
                                                ])
                # print("S vale: ", S)
                # exit()
                for ino in range(npe):
                    for idir in range(ngl):
                        #cálculo do Fint(i)
                        dA1dyi = np.array([
                                           [dfidksi[ino]*(1-idir),dfideta[ino]*(1-idir)],
                                           [dfidksi[ino]*(idir),dfideta[ino]*(idir)]  
                                          ])
                        # print("DA_ab vale: ", dA1dyi)
                        
                        aux = np.dot(A1,np.linalg.inv(A0))
                        aux = np.dot(np.transpose(dA1dyi),aux)
                        aux1 = np.dot(np.transpose(np.linalg.inv(A0)),aux)
                        # print("DE_1 vale: ", aux1)
                        
                        aux = np.dot(dA1dyi,np.linalg.inv(A0))
                        aux = np.dot(np.transpose(A1),aux)
                        aux2 = np.dot(np.transpose(np.linalg.inv(A0)),aux)
                        # print("DE_2 vale: ", aux2)
                        dEdyi = 1/2 * (aux1 + aux2)
                        # print("DE_ab vale: ", dEdyi)
                        
                        ef[ngl * ino + idir] = h * np.tensordot(S,dEdyi,axes=2)
                        print("f_int_i vale: ", ef)
                        
                        for jno in range(npe):
                            for jdir in range(ngl):
                                #cálculo H(i,j)    
                                dA1dyj = np.array([
                                                   [dfidksi[jno]*(1-jdir),dfideta[jno]*(1-jdir)],
                                                   [dfidksi[jno]*(jdir),dfideta[jno]*(jdir)]  
                                                  ])

                                # print("DA1_gz vale: ", dA1dyj)
                                
                                aux = np.dot(A1,np.linalg.inv(A0))
                                aux = np.dot(np.transpose(dA1dyj),aux)
                                aux1 = np.dot(np.transpose(np.linalg.inv(A0)),aux)

                                aux = np.dot(dA1dyj,np.linalg.inv(A0))
                                aux = np.dot(np.transpose(A1),aux)
                                aux2 = np.dot(np.transpose(np.linalg.inv(A0)),aux)

                                dEdyj = 1/2 * (aux1 + aux2)
                                # counter+=1 
                                # print("Counter:",counter)
                                # print("DE_gz vale: ", dEdyj)
                                # print("_______________________")
                                
                                dSdEdEdyj = young / (1 - ni**2) * np.array([
                                                  [dEdyj[0,0] + ni * dEdyj[1,1],        (1 - ni) * dEdyj[0,1]],
                                                  [       (1 - ni) * dEdyj[1,0], ni * dEdyj[0,0] + dEdyj[1,1]]
                                                ])
                                
                                # counter+=1 
                                # print("Counter:",counter)
                                # # print("DS_gz vale: ", dSdEdEdyj)
                                # print("_______________________")
                                aux = np.dot(dA1dyj,np.linalg.inv(A0))
                                aux = np.dot(np.transpose(dA1dyi),aux)
                                aux1 = np.dot(np.transpose(np.linalg.inv(A0)),aux)

                                aux = np.dot(dA1dyi,np.linalg.inv(A0))
                                aux = np.dot(np.transpose(dA1dyj),aux)
                                aux2 = np.dot(np.transpose(np.linalg.inv(A0)),aux)

                                d2Edyidyj = 1/2 * (aux1 + aux2)
                                # print('D2E',d2Edyidyj)
                                # exit()
                                h_abgz = 0.0
                                for i in range(2):
                                    for j in range(2):
                                        h_abgz=h_abgz + h *( dSdEdEdyj[i,j] * dEdyi[i,j] + S[i,j] * d2Edyidyj[i,j])
                                        # print(f'i: {i} j: {j} DE_gz: {dEdyi[i,j]} DS_gz: {d2Edyidyj[i,j]} S_Piola: {S[i,j]} D2E: {dSdEdEdyj[i,j]}')
                                # print(f"h_abgz:  {h_abgz}")
                                # print("iel:",iel)
                                # print("idir:",idir)
                                # print("jdir:",jdir)
                                # print("ino:",ino)
                                # print("jno:",jno)
                                # print(f"i:{ngl * conec[iel,ino] + idir}     j:{ngl * conec[iel,jno] + jdir}")
                                # print(f"ELEMS(iel,ino):{conec[iel,ino]}")
                                # print(f"ELEMS(iel,jno):{conec[iel,jno]}")
                                # exit()
                                # print('eh',h * (np.tensordot(dSdEdEdyj,dEdyi) + np.tensordot(S,d2Edyidyj)))
                                # print("eh:")
                                # print("_________")
                                # if ngl * ino + idir != ngl * jno + jdir:
                                #     print('eh out of diagonal',h * (np.tensordot(dSdEdEdyj,dEdyi) + np.tensordot(S,d2Edyidyj)))

                                eh[ngl * ino + idir,ngl * jno + jdir] = h_abgz
                                # eh[ngl * ino + idir,ngl * jno + jdir] = h * (np.tensordot(dSdEdEdyj,dEdyi) + np.tensordot(S,d2Edyidyj))


                                
                for ino in range(npe):
                    for idir in range(ngl):
                        fint[ngl * conec[iel,ino] + idir] += ef[ngl * ino + idir] * peso * jota
                        fext[ngl * conec[iel,ino] + idir] += efext[ngl * ino + idir] * peso * jota
                        for jno in range(npe):
                            for jdir in range(ngl):
                                hessiana[ngl * conec[iel,ino] + idir,ngl * conec[iel,jno] + jdir] += eh[ngl * ino + idir,ngl * jno + jdir] * peso * jota
                

            # print("mfifi vale:", mfifi)
            # print("mfi vale:", mfi)
            # print("A0 vale:", A0)
            # exit()
        #contribuição das fibras
        for ifibras in range(nfibras):
            young = fproperties[fconec[ifibras,2],0]
            area = fproperties[fconec[ifibras,2],1]
            x11 = fnodes[fconec[ifibras,0],0]
            x21 = fnodes[fconec[ifibras,0],1]
            x12 = fnodes[fconec[ifibras,1],0]
            x22 = fnodes[fconec[ifibras,1],1]
            L0 = ((x12-x11)**2.0+(x22-x21)**2.0)**0.5
            y11 = yf[fconec[ifibras,0],0]
            y21 = yf[fconec[ifibras,0],1]
            y12 = yf[fconec[ifibras,1],0]
            y22 = yf[fconec[ifibras,1],1]
            l = ((y12-y11)**2.0+(y22-y21)**2.0)**0.5
            
            k1i = ifconec[fconec[ifibras,0]][1]
            k2i = ifconec[fconec[ifibras,0]][2]
            k1f = ifconec[fconec[ifibras,1]][1]
            k2f = ifconec[fconec[ifibras,1]][2]
            # print("ifconec",ifconec)
            # print("k1i, k2i, k1f, k2f",k1i, k2i, k1f, k2f )
            # exit()

            vfi,vdfdk1i,vdfdk2i = fforma(k1i,k2i)
            vff,vdfdk1f,vdfdk2f = fforma(k1f,k2f)

            # print ("vfi, vff", vfi, vff)
            # exit()


            mphi = np.zeros((4,2*ngl*npe))
            for j in range(npe):
                for k in range(ngl):
                    mphi[0,ngl*j+k] = vfi[j]*(1-k)
                    mphi[1,ngl*j+k] = vfi[j]*(k)
                    mphi[2,ngl*npe+ngl*j+k] = vff[j]*(1-k)
                    mphi[3,ngl*npe+ngl*j+k] = vff[j]*(k)
            # print("mphi vale: ", mphi)
            # exit()

            ifelem = [ifconec[fconec[ifibras,0]][0],ifconec[fconec[ifibras,1]][0]]
            index = []
            # print("Ifelem vale: ", ifelem)

            for i in range(2):
                for j in range(npe):
                    for k in range(ngl):
                        index.append(ngl*conec[ifelem[i],j]+k)               
            # print("index vale: ", index)
            # exit()
            Eg = 1.0/2.0 * (l**2.0 - L0**2.0) / L0**2.0
            S = young * Eg
            eff = np.array([
                           -((area*S*(-y11 + y12))/L0),
                           -((area*S*(-y21 + y22))/L0),
                           (area*S*(-y11 + y12))/L0,
                           (area*S*(-y21 + y22))/L0
                           ])
            ehf = np.array([
                            [
                             (area*S)/L0 + (area*(y11 - y12)**2* young)/L0**3,
                             (area*(y11 - y12)* (y21 - y22)*young)/L0**3,
                             -((area*S)/L0) - (area*(y11 - y12)**2*young)/L0**3,
                             -((area*(y11 - y12)*(y21 - y22)*young)/L0**3)
                            ],
                            [
                             (area*(y11 - y12)*(y21 - y22)*young)/L0**3,
                             (area*S)/L0 + (area*(y21 - y22)**2*young)/L0**3,
                             -((area*(y11 - y12)*(y21 - y22)*young)/L0**3),
                             -((area*S)/L0) - (area*(y21 - y22)**2*young)/L0**3
                            ],
                            [
                             -((area*S)/L0) - (area*(y11 - y12)**2*young)/L0**3,
                             -((area*(y11 - y12)*(y21 - y22)*young)/L0**3),
                             (area*S)/L0 + (area*(y11 - y12)**2*young)/L0**3,
                             (area*(y11 - y12)*(y21 - y22)*young)/L0**3
                            ],
                            [
                             -((area*(y11 - y12)*(y21 - y22)*young)/L0**3),
                             -((area*S)/L0) - (area*(y21 - y22)**2*young)/L0**3,
                             (area*(y11 - y12)*(y21 - y22)*young)/L0**3,
                             (area*S)/L0 + (area*(y21 - y22)**2*young)/L0**3
                            ]
                           ])
            ff = np.dot(np.transpose(mphi),eff)
            hf = np.dot(np.dot(np.transpose(mphi),ehf),mphi)

            # print("ff vale: ", ff)
            # print("hf vale: ", hf)
            # exit()

            for i in range(2*ngl*npe):
                fint[index[i]] += ff[i]
                for j in range(2*ngl*npe):
                    hessiana[index[i],index[j]] += hf[i,j]
            
     
        # print("hessiana vale: ", hessiana)
        # print("fint vale: ", fint)
        # exit()
        g = - fint + fext

        # print("g vale: ", g)
        # exit()

        #aplicação das condições de contorno
        for ivinc in vinc:
            node = int(ivinc[0])
            dir = int(ivinc[1])
            val = ivinc[2]
            g[ngl * node + dir] = 0.0
            hessiana[ngl * node + dir,:] = 0.0
            hessiana[:,ngl * node + dir] = 0.0
            hessiana[ngl * node + dir,ngl * node + dir] = 1.0


        # print("g vale: ", g)
        # print("hessiana vale: ", hessiana)
        # exit()
        #resolver o sistema
        
        delta_Y = np.linalg.solve(hessiana,g)


        #atualizar as posições finais dos nós
        for ino in range(nnodes):
            for igl in range(ngl):
                y[ino,igl] += delta_Y[ngl * ino + igl]

        # print("delta_y atualizado vale: ", delta_Y)
        # print("y atualizado vale: ", y)
        # exit()


        #atualizar as posições finais dos nós das fibras
        yf *= 0.0
        for ifn in range(nnodesf):
            ielem = ifconec[ifn][0]
            ksi1 = ifconec[ifn][1]
            ksi2 = ifconec[ifn][2]
            vf,vdfdk1i,vdfdk2i = fforma(ksi1,ksi2)
            
            for i in range(npe):
                for j in range(ngl):
                    yf[ifn,j] += vf[i] * y[conec[ielem,i],j]

        # print("yf ", yf)
        # exit()
        erro = np.linalg.norm(delta_Y) / np.linalg.norm(nodes.reshape(-1,1))
        if iter == 1:
            print("hesiana vale: ", hessiana)
            print("erro: ",erro)
            print("fint vale: ", fint)
            print("fext vale: ", fext)
            print("g vale: ", g)
            print("delta_Y vale: ", delta_Y)
            # exit()
        iter += 1
        print("     iter: ",iter," | erro: ",f'{erro:.08f}')
        print("delta_Y vale: ", delta_Y)


    #cálculo das tensões
    sigma = np.zeros((nnodes,3))
    for iel in range(nelems):
        young = properties[conec[iel,-1],0]
        ni = properties[conec[iel,-1],1]
        h = properties[conec[iel,-1],2]

        for ih in range(npe):
            ksi = ksieta[ih,0]
            eta = ksieta[ih,1]
            fi,dfidksi,dfideta = fforma(ksi,eta)
            A0 = np.zeros((2,2)) 
            A1 = np.zeros((2,2)) 
            for ino in range(npe):
                A0[0,0] += dfidksi[ino] * nodes[conec[iel,ino],0]
                A0[0,1] += dfideta[ino] * nodes[conec[iel,ino],0]
                A0[1,0] += dfidksi[ino] * nodes[conec[iel,ino],1]
                A0[1,1] += dfideta[ino] * nodes[conec[iel,ino],1]
                A1[0,0] += dfidksi[ino] * y[conec[iel,ino],0]
                A1[0,1] += dfideta[ino] * y[conec[iel,ino],0]
                A1[1,0] += dfidksi[ino] * y[conec[iel,ino],1]
                A1[1,1] += dfideta[ino] * y[conec[iel,ino],1]
            
            jota = np.linalg.det(A0)

            #EPT
            A = np.dot(A1,np.linalg.inv(A0))
            C = np.dot(np.transpose(A),A)
            E = 1/2 * (C - np.identity(2))
            S = young / (1 - ni**2) * np.array([
                                                [E[0,0] + ni * E[1,1],    (1 - ni) * E[0,1]],
                                                [   (1 - ni) * E[1,0], ni * E[0,0] + E[1,1]]
                                            ])
            
            sig = S#1.0 / jota / h * np.dot(A,np.dot(S,np.transpose(A)))

            sigma[conec[iel,ih],0] += sig[0,0] / iconec[conec[iel,ih]]
            sigma[conec[iel,ih],1] += sig[1,1] / iconec[conec[iel,ih]]
            sigma[conec[iel,ih],2] += sig[0,1] / iconec[conec[iel,ih]]

    saida.append('#\n')
    saida.append('desloc.x\n')
    for i in range(nnodes):
        saida.append(f'{y[i,0]-nodes[i,0]} {y[i,1]-nodes[i,1]} 0 {y[i,0]-nodes[i,0]}\n')
    saida.append('#\n')
    saida.append('desloc.y\n')
    for i in range(nnodes):
        saida.append(f'{y[i,0]-nodes[i,0]} {y[i,1]-nodes[i,1]} 0 {y[i,1]-nodes[i,1]}\n')
    saida.append('#\n')
    saida.append('sig.x\n')
    for i in range(nnodes):
        saida.append(f'{y[i,0]-nodes[i,0]} {y[i,1]-nodes[i,1]} 0 {sigma[i,0]}\n')
    saida.append('#\n')
    saida.append('sig.y\n')
    for i in range(nnodes):
        saida.append(f'{y[i,0]-nodes[i,0]} {y[i,1]-nodes[i,1]} 0 {sigma[i,1]}\n')
    saida.append('#\n')
    saida.append('tal.xy\n')
    for i in range(nnodes):
        saida.append(f'{y[i,0]-nodes[i,0]} {y[i,1]-nodes[i,1]} 0 {sigma[i,2]}\n')

arquivo = open('saida.ogl','w')
arquivo.writelines(saida)
arquivo.close()
print("fim da análise")