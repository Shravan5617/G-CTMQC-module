  )-  ^   k820309    ?          19.0        Սbb                                                                                                          
       time_evolution.f90 TIME_EVOLUTION          @       ?                                  
                        ?                                  
                                                           
                                                           
                                                           
                                                           
                                                           
                                                           
                                                      	     
                                                    
     KIND                                                                                                                                                                                                             #         @                                                      #BOSIGMA    #RCL    #VCL    #TIME             
                                                             p        5 r    p        5 r    p          5 r      5 r      5 r        5 r      5 r      5 r                               
                                                     
      p        5 r    p          5 r      5 r        5 r      5 r                               
                                                     
      p        5 r    p          5 r      5 r        5 r      5 r                                
                                                       D                                             #         @                                                                  @                                                                                                               D  @                                                     D                                      
                 D                                                        D                                             #         @                                                      #COEFF    #FORCE    #TRAJLABEL             
                                                         p          5 r        5 r                               
                                                    
       p        5 r    p          5 r      5 r        5 r      5 r                                
                                             #         @                                                       #V !   #FORCE "            
                                !                    
     p          5 r        5 r                               
                                 "                    
    p          5 r        5 r                                                                #                   
                &                   &                   &                                                                                      $                                   &                                           #         @                                  %                    #X &   #V '            
                                &                    
     p          5 r        5 r                               
                                 '                    
    p          5 r        5 r                      #         @                                  (                    #Q )   #TRAJLABEL *            
                                 )                    
    p          5 r        5 r                                
                                  *           #         @                                  +                    #V ,   #COEFF -   #K_LI .   #TRAJLABEL /            
                                 ,                    
    p          5 r        5 r                               
                                -                         p          5 r        5 r                               
                                 .                    
      p        5 r    p          5 r      5 r        5 r      5 r                                
                                  /           #         @                                  0                    #COEFF 1   #V 2   #TRAJLABEL 3            
                                1                     
    p          5 r        5 r                               
                                 2                    
    p          5 r        5 r                                                                 3            #         @                                  4                    #MY_RHO 5   #V 6   #R 7   #TRAJLABEL 8            
                                 5                          p        5 r    p          5 r      5 r        5 r      5 r                               
                                6                    
     p          5 r        5 r                               
                                 7                    
    p          5 r        5 r                                
                                  8           #         @                                  9                    #COEFF :   #FORCE ;   #ACC_FORCE <   #K_LI =   #TRAJLABEL >            
                                 :                        p          5 r        5 r                               
                                ;                    
 
    p          5 r        5 r                               
                                 <                    
      p        5 r    p          5 r      5 r        5 r      5 r                               
                                 =                    
 	     p        5 r    p          5 r      5 r        5 r      5 r                                
                                  >                     D                                              #         @                                  ?                    #RCL @   #ACC_FORCE C   #BOSIGMA E   #K_LI F            
                                 @                    
      p        5 r A   p          5 r A     5 r B       5 r A     5 r B                              
                                 C                    
        p        5 r B   p        5 r A   p          5 r A     5 r B     5 r D       5 r A     5 r B     5 r D                              
                                 E                            p        5 r D   p        5 r A   p          5 r A     5 r D     5 r D       5 r A     5 r D     5 r D                              
                                F                    
         p        5 r D   p        5 r A   p          5 r A     5 r D     5 r D       5 r A     5 r D     5 r D                               D  @                              G     d                 D                                 H     
                                                  I                                   &                                                                                     J                   
                &                                                                                     K                   
                &                                                     D                                 L     
                 D                                                        D                                  M                                                      N                   
                &                                                                                     O                   
                &                   &                                                                                     P                   
                &                   &                                                                                     Q                   
                &                                                    @ @                              R                   
                &                   &                                                    @ @                              S                   
                &                   &                                                    @ @                              T                   
                &                                                    @ @                              U                   
                &                   &                   &                                                    @ @                              V                   
                &                   &                   &                                                                                     W                   
                &                   &                                                    @ @                              X                                   &                   &                   &                                                    @ @                              Y                                   &                   &                                           #         @                                   Z                     #         @                                  [                     #         @                                  \                     #         @                                  ]                               D                                 D                      D                                 A                      D                                 B               ?   *      fn#fn    ?   @   J   VARIABLES    
  @   J   KINDS #   J  @   J   ELECTRONIC_PROBLEM '   ?  @   J   COEFFICIENTS_EVOLUTION $   ?  @   J   CLASSICAL_EVOLUTION &   
  @   J   COHERENCE_CORRECTIONS    J  @   J   OUTPUT    ?  @   J   SHOPPING '   ?  @   J   TRAJECTORIES_SELECTION    
  =       KIND+KINDS    G  p       DP+KINDS    ?  p       QP+KINDS    '  q       PLOT+OUTPUT $   ?    a   PLOT%BOSIGMA+OUTPUT     ?  ?   a   PLOT%RCL+OUTPUT     ?  ?   a   PLOT%VCL+OUTPUT !   T  @   a   PLOT%TIME+OUTPUT "   ?  @       TYP_CAL+VARIABLES /   ?  H       GENERATE_RANDOM_SEED_HOP+TOOLS %     @       COUNT_TRAJ+VARIABLES !   \  @       NSTEPS+VARIABLES    ?  @       DUMP+VARIABLES    ?  @       DT+VARIABLES     	  @       NTRAJ+VARIABLES $   \	  @       TYPE_DECO+VARIABLES :   ?	  m       ACCUMULATED_BOFORCE+COHERENCE_CORRECTIONS @   	
  ?   a   ACCUMULATED_BOFORCE%COEFF+COHERENCE_CORRECTIONS @   ?
  ?   a   ACCUMULATED_BOFORCE%FORCE+COHERENCE_CORRECTIONS D   q  @   a   ACCUMULATED_BOFORCE%TRAJLABEL+COHERENCE_CORRECTIONS 4   ?  Z       UPDATE_VELOCITY+CLASSICAL_EVOLUTION 6     ?   a   UPDATE_VELOCITY%V+CLASSICAL_EVOLUTION :   ?  ?   a   UPDATE_VELOCITY%FORCE+CLASSICAL_EVOLUTION "   3  ?       BOFORCE+VARIABLES $   ?  ?       OCC_STATE+VARIABLES 4   {  V       UPDATE_POSITION+CLASSICAL_EVOLUTION 6   ?  ?   a   UPDATE_POSITION%X+CLASSICAL_EVOLUTION 6   e  ?   a   UPDATE_POSITION%V+CLASSICAL_EVOLUTION -   ?  ^       BOPROBLEM+ELECTRONIC_PROBLEM /   W  ?   a   BOPROBLEM%Q+ELECTRONIC_PROBLEM 7   ?  @   a   BOPROBLEM%TRAJLABEL+ELECTRONIC_PROBLEM 1   +  s       RK4_COEFF+COEFFICIENTS_EVOLUTION 3   ?  ?   a   RK4_COEFF%V+COEFFICIENTS_EVOLUTION 7   2  ?   a   RK4_COEFF%COEFF+COEFFICIENTS_EVOLUTION 6   ?  ?   a   RK4_COEFF%K_LI+COEFFICIENTS_EVOLUTION ;   ?  @   a   RK4_COEFF%TRAJLABEL+COEFFICIENTS_EVOLUTION 0   ?  i       DECOHERENCE_COORECTION+SHOPPING 6   C  ?   a   DECOHERENCE_COORECTION%COEFF+SHOPPING 2   ?  ?   a   DECOHERENCE_COORECTION%V+SHOPPING :   k  @   a   DECOHERENCE_COORECTION%TRAJLABEL+SHOPPING !   ?  q       HOPPING+SHOPPING (     ?   a   HOPPING%MY_RHO+SHOPPING #   ?  ?   a   HOPPING%V+SHOPPING #   ?  ?   a   HOPPING%R+SHOPPING +     @   a   HOPPING%TRAJLABEL+SHOPPING 8   X  ?       NON_ADIABATIC_FORCE+CLASSICAL_EVOLUTION >   ?  ?   a   NON_ADIABATIC_FORCE%COEFF+CLASSICAL_EVOLUTION >   r  ?   a   NON_ADIABATIC_FORCE%FORCE+CLASSICAL_EVOLUTION B     ?   a   NON_ADIABATIC_FORCE%ACC_FORCE+CLASSICAL_EVOLUTION =   ?  ?   a   NON_ADIABATIC_FORCE%K_LI+CLASSICAL_EVOLUTION B   ?  @   a   NON_ADIABATIC_FORCE%TRAJLABEL+CLASSICAL_EVOLUTION "   ?  @       NSTATES+VARIABLES 7   .  w       QUANTUM_MOMENTUM+COHERENCE_CORRECTIONS ;   ?  ?   a   QUANTUM_MOMENTUM%RCL+COHERENCE_CORRECTIONS A   y    a   QUANTUM_MOMENTUM%ACC_FORCE+COHERENCE_CORRECTIONS ?   ?    a   QUANTUM_MOMENTUM%BOSIGMA+COHERENCE_CORRECTIONS <   ?    a   QUANTUM_MOMENTUM%K_LI+COHERENCE_CORRECTIONS *   ?   @       MODEL_POTENTIAL+VARIABLES &   ?   @       C_PARAMETER+VARIABLES *   5!  ?       INITIAL_BOSTATE+VARIABLES    ?!  ?       R0+VARIABLES    M"  ?       K0+VARIABLES %   ?"  @       FINAL_TIME+VARIABLES     #  @       N_DOF+VARIABLES $   Y#  @       N_INIT_BO+VARIABLES )   ?#  ?       WEIGHT_BOSTATE+VARIABLES ,   %$  ?       INITIAL_POSITIONS+VARIABLES *   ?$  ?       INITIAL_MOMENTA+VARIABLES    m%  ?       MASS+VARIABLES    ?%  ?       RCL    ?&  ?       VCL     A'  ?       CLASSICAL_FORCE    ?'  ?       MY_FORCE    ?(  ?       K_LI    E)  ?       TDVP    ?)  ?       BOSIGMA    ?*  ?       BOCOEFF    I+  H       EVOLUTION    ?+  H       INPUT_SUMMARY &   ?+  H       INITIALIZE_LOCAL_VARS $   !,  H       FINALIZE_LOCAL_VARS 3   i,  @     QUANTUM_MOMENTUM%NSTATES+VARIABLES 1   ?,  @     QUANTUM_MOMENTUM%NTRAJ+VARIABLES 1   ?,  @     QUANTUM_MOMENTUM%N_DOF+VARIABLES 