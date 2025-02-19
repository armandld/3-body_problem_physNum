#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template
#include <numeric>

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{
private:

// EngineEuler specific members
  unsigned int maxit; // nombre maximale d iterations
  double tol;         // tolerance methode iterative
  double alpha;       // parametre pour le scheme d'Euler

// définition des variables
double tfin;         // Temps final
unsigned int nsteps; // Nombre de pas de temps
double ml;           // Masse de la Lune
double mt;           // Masse de la Terre
double dist;         // Distance Terre-Lune
double Om;           // Vitesse de rotation du repère
double G_grav;       // Constante gravitationnelle
double xt;           // Position de la Terre
double xl;           // Position de la Lune
double dist_s_t;     // Distance satellite-Terre
double dist_s_l;     // Distance satellite-Lune

  valarray<double> y0 = std::valarray<double>(0.e0, 4); // Correctly initialized
  valarray<double> y  = std::valarray<double>(0.e0, 4); // Correctly initialized

  double t,dt;  // Temps courant pas de temps

  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  void printOut(bool write)
  {
  // TODO calculer l'energie mecanique
    double Energy =  1/2*(y[0]*y[0]+y[1]*y[1])+G_grav*mt/dist_s_t+G_grav*ml/dist_s_l;

    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      *outputFile << t << " " << y[0] << " " << y[1] << " " \
      << y[2] << " " << y[3] << " " << Energy << " "<< endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

    void compute_f(valarray<double>& f) //  TODO: Calcule le tableau de fonctions f(y)
    {

		
		f[2] = f[0]; // vitesse en x du satellite
		f[3] = f[1]; // vitesse en y du satellite
		
		f[0] = pow(Om, 2) * f[2] + 2 * Om * f[1] 
       - G_grav * mt / pow(dist_s_t,2) * f[2] / sqrt(pow(f[2], 2) + pow(f[3], 2)) 
       - G_grav * ml / pow(dist_s_l,2) * f[2] / sqrt(pow(f[2], 2) + pow(f[3], 2));

		f[3] = pow(Om, 2) * f[1] - 2 * Om * f[0] 
       - G_grav * mt / pow(dist_s_t, 2) * f[3] / sqrt(pow(f[2], 2) + pow(f[3], 2)) 
       - G_grav * ml / pow(dist_s_l, 2) * f[3] / sqrt(pow(f[2], 2) + pow(f[3], 2));


    }

    // New step method from EngineEuler
    void step()
    {
      unsigned int iteration=0;
      double error=999e0;
      valarray<double> f =valarray<double>(0.e0,4); 
      valarray<double> yold=valarray<double>(y);
      valarray<double> y_control=valarray<double>(y);
      valarray<double> y_previous=valarray<double>(y);
      valarray<double> delta_y_EE=valarray<double>(y);

	  dist_s_l = sqrt(pow(y[2]-xl,2)+pow(y[3],2));
	  dist_s_t = sqrt(pow(y[2]-xt,2)+pow(y[3],2));

      //TODO : écrire un algorithme valide pour chaque alpha dans [0,1]
      // tel que alpha=1 correspond à Euler explicite et alpha=0 à Euler implicite 
      // et alpha=0.5 à Euler semi-implicite
      if(alpha >= 0. && alpha <= 1.0){
        t += dt;                 //mise à jour du temps 
        valarray<double> yold_copie=valarray<double>(yold); //copie de yold
        compute_f(yold); 				  //yold vaut maintenant sa dérivée
        
        while(error>tol && iteration<=maxit)
        {
			valarray<double> y_copie=valarray<double>(y); //copie de y
			compute_f(y);				//y vaut maintenant sa dérivée
            y = yold_copie + (alpha * yold + (1-alpha) * y) * dt;
            valarray<double> y_copie_bis=valarray<double>(y);
            compute_f(y);
            error = abs((y_copie_bis - yold_copie - (alpha * yold + (1 - alpha)) * y ) * dt).max();
            y = y_copie_bis;
            ++iteration;
		}	
        if(iteration>=maxit){
          cout << "WARNING: maximum number of iterations reached, error: " << error << endl;
        }
      }
      else
      {
        cerr << "alpha not valid" << endl;
      }
      cout << iteration << endl;
  
    }

public:
    // Modified constructor
    Engine(ConfigFile configFile)
    {
      // Stockage des parametres de simulation dans les attributs de la classe
      tfin     = configFile.get<double>("tfin",tfin);	        // lire le temps final de simulation
      nsteps   = configFile.get<unsigned int>("nsteps",nsteps); // lire le nombre de pas de temps
      y0[0]    = configFile.get<double>("vx0",y0[0]);  // vitesse initiale selon x	    
      y0[1]    = configFile.get<double>("vy0",y0[1]);  // vitesse initiale selon y       
      y0[2]    = configFile.get<double>("x0",y0[2]);   // position initiale selon x       
      y0[3]    = configFile.get<double>("y0",y0[3]);   // position initiale selon y	    
      G_grav   = configFile.get<double>("G_grav",G_grav);           
      ml       = configFile.get<double>("ml",ml);            
      mt       = configFile.get<double>("mt",mt);        
      dist     = configFile.get<double>("dist",dist);        
      sampling = configFile.get<unsigned int>("sampling",sampling);
      tol      = configFile.get<double>("tol", tol);
      maxit    = configFile.get<unsigned int>("maxit", maxit);
      alpha    = configFile.get<double>("alpha", alpha);
      // TODO: calculer le time step
      dt       = tfin/nsteps;

      
      // Ouverture du fichier de sortie
      outputFile = new ofstream(configFile.get<string>("output","output.out").c_str()); 
      outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
    };


    // Destructeur virtuel
    virtual ~Engine()
    {
      outputFile->close();
      delete outputFile;
    };
      // Simulation complete
    void run()
    {
      // TODO : initialiser la position de la Terre et de la Lune, ainsi que la position de X' du satellite et Omega
      Om = sqrt(G_grav*mt/(xl*dist*dist));
      xt = -ml*dist/(mt+ml);
      xl = mt*dist/(mt+ml); 
      y0[2] = (xt*sqrt(ml)+xl*sqrt(mt))/(sqrt(mt)+sqrt(ml));
      

      t = 0.e0; // initialiser le temps
      y = y0;   // initialiser la position 
      last = 0; // initialise le parametre d'ecriture

      printOut(true); // ecrire la condition initiale

      for(unsigned int i(0); i<nsteps; ++i) // boucle sur les pas de temps
      {
        step();  // faire un pas de temps
        printOut(false); // ecrire le pas de temps actuel
      }
      printOut(true); // ecrire le dernier pas de temps

    };
   
};

// programme
int main(int argc, char* argv[])
{
  string inputPath("configuration.in.example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

  Engine* engine;

  // Create an instance of Engine instead of EngineEuler
  engine = new Engine(configFile);

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}




