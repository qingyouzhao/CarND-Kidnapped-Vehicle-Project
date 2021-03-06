/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <assert.h>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // todo: I need a better number?
  constexpr int my_mum_particles = 1000;
  num_particles = my_mum_particles;  // TODO: Set the number of particles
  particles.clear();

  std::default_random_engine gen;
  const double& std_x     = std[0];
  const double& std_y     = std[1];
  const double& std_theta = std[2];
  std::normal_distribution<double> dist_x(x, std_x);
  std::normal_distribution<double> dist_y(y, std_y);
  std::normal_distribution<double> dist_theta(theta, std_theta);

  for (unsigned int i =0; i < num_particles; i++)
  {
    Particle p;
    p.id      = i;
    p.x       = dist_x(gen);
    p.y       = dist_y(gen);
    p.theta   = dist_theta(gen);
    p.weight  = 1.0;
    // PrintParticle(p);
    // I can confirm particles are initialized fine.
    // todo: do I need to init the associations here?
    particles.push_back(p);
  }
  assert(particles.size() == num_particles);
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[],
  double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   // prepare gen
  std::default_random_engine gen;
  const double& std_x = std_pos[0];
  const double& std_y = std_pos[1];
  const double& std_theta = std_pos[2];

  //std::cout << "Prediction step" << std::endl;
  //std::cout << "delta_t = " << delta_t << std::endl;
  //std::cout << "velocity = " << velocity << std::endl;
  //std::cout << "yaw_rate = " << yaw_rate << std::endl;


  //std::cout << "std_x = "<< std_x  << std::endl;
  //std::cout << "std_y = " << std_y << std::endl;
  //std::cout << "std_theta = "<< std_theta << std::endl;


  std::normal_distribution<double> dist_x(0.0, std_x);
  std::normal_distribution<double> dist_y(0.0, std_y);
  std::normal_distribution<double> dist_theta(0.0, std_theta);

  // Predict new state based on v and yaw_rate;
  for (Particle& p : particles)
  {
    assert(!std::isnan(p.x));
    assert(!std::isnan(p.y));
    //std::cout << "before prediction step" << std::endl;
    //PrintParticle(p);
    // Retrieve old states
    const double x0 = p.x;
    const double y0 = p.y;
    const double theta0 = p.theta;
    // Calculate new states
    const double delta_yaw = yaw_rate * delta_t;
    const double thetaf = theta0 + delta_yaw;
    const double xf = x0 + velocity / yaw_rate * (sin(thetaf) - sin(theta0));
    const double yf = y0 + velocity / yaw_rate * (cos(theta0) - cos(thetaf));
    // Update new states to particle and add noise
    // todo: verify this is the correct noise we want to add.
    p.x = xf + dist_x(gen);
    p.y = yf + dist_y(gen);
    if( (std::isnan(p.x)) || (std::isnan(p.y)))
    {
      p.x = x0;
      p.y = y0;
    }
    p.theta = thetaf + dist_theta(gen);
    //std::cout << "After prediction step" << std::endl;
    //PrintParticle(p);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  // Match observations with predicted?
  // Find predicted measurement. 
  for (LandmarkObs& ob : observations)
  {
    // find closest observation
    double min_dist_sqr = std::numeric_limits<double>::max();
    for (const LandmarkObs& prediction : predicted)
    {
      auto GetDistSqr = [](LandmarkObs ob1, LandmarkObs ob2) {return (ob1.x - ob2.x) *(ob1.x - ob2.x) + (ob1.y - ob2.y)*(ob1.y - ob2.y); };
      const double dist_sqr = GetDistSqr(prediction, ob);
      if (min_dist_sqr > dist_sqr)
      {
        min_dist_sqr = dist_sqr;
        ob.id = prediction.id;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  // todo: look at update weight algo again.
  // 
  // todo: What does sensor range has to do with our observation.
  // todo: What is std_landmark x, y in? map space?

  // First we construct predicted LandmarkObs from landmarks in map space
  // std::cout << "Begin " << "updateWeights" << std::endl;
  assert(observations.size() != 0);
  vector<LandmarkObs> predicted;
  for (const Map::single_landmark_s lm : map_landmarks.landmark_list)
  {
    LandmarkObs predicted_lm;
    predicted_lm.id = lm.id_i;
    predicted_lm.x = lm.x_f;
    predicted_lm.y = lm.y_f;
    // PrintLandmarkObs(predicted_lm);
    predicted.emplace_back(predicted_lm);
  }

  // Then we update weight for each particle
  for(Particle& p : particles)
  {
    assert(!std::isnan(p.x));
    assert(!std::isnan(p.y));

    // p is in Map coordinate system.
    // Step 1: convert observations to map spaces based on p position
    vector<LandmarkObs> observation_ms;
    for (const LandmarkObs& obs_ps : observations)
    {
      LandmarkObs obs_ms = GetMapSpaceObservation(p, obs_ps);
      observation_ms.push_back(obs_ms);
    }

    
    // Step2: Match observations with predicted, the observation id will be updated
    dataAssociation(predicted, observation_ms);

    // Step3: update weights based on observation and prediction
    double final_weight = 1.0;
    for (const LandmarkObs& obs : observation_ms)
    {
      auto IsMatchingObs = [&obs](LandmarkObs pred) {return pred.id == obs.id; };
      LandmarkObs matching_prediction = *std::find_if(predicted.begin(), predicted.end(), IsMatchingObs);
      double sig_x = std_landmark[0];
      double sig_y = std_landmark[1];
      double x_obs = obs.x; // Map space of obs
      double y_obs = obs.y;
      double mu_x = matching_prediction.x; // map space of LM
      double mu_y = matching_prediction.y; 
      double prob = multiv_prob(sig_x,sig_y,x_obs,y_obs,mu_x,mu_y);
      if (prob < 0.0)
      {
        PrintLandmarkObs(obs);
        PrintParticle(p);
        assert(false);
      }
      final_weight *= prob;
    }
    p.weight = final_weight;
    assert(p.weight >= 0.0);
    // Then we have to update weights
  }  
  // /std::cout << "End " << "updateWeights" << std::endl;
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  // std::cout << "Begin " << "resample" << std::endl;

  weights.clear();
  assert(particles.size() == num_particles);
  for (unsigned int i = 0; i < particles.size(); i++)
  {
    weights.push_back(particles[i].weight);
  }
  std::discrete_distribution<int> weight_dist(weights.begin(), weights.end());
  std::default_random_engine gen;

  // Set of current particles
  std::vector<Particle> resampled_particles;
  for (unsigned int n = 0; n < num_particles; n++)
  {
    int sample_index = weight_dist(gen);
    assert(particles[sample_index].weight >= 0.0);
    resampled_particles.push_back(particles[sample_index]);
  }

  assert(resampled_particles.size() == num_particles);
  particles.swap(resampled_particles);
  assert(particles.size() == num_particles);
  // std::cout << "End " << "resample" << std::endl;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

LandmarkObs ParticleFilter::GetMapSpaceObservation(const Particle& p, const LandmarkObs& obs_ps)
{
  /*
      x        y
      |        |
  y___|        |____x
      Car     Map
*/

  LandmarkObs obs_ms; // map space wrt to current particle
  obs_ms.x = p.x + cos(p.theta) * obs_ps.x - sin(p.theta) * obs_ps.y;
  obs_ms.y = p.y + sin(p.theta) * obs_ps.x + cos(p.theta) * obs_ps.y;
  return obs_ms;
}

void ParticleFilter::PrintParticle(const Particle& p) const
{
  std::cout << "Particle: { " << std::endl;
  std::cout << "id: " << p.id << std::endl;
  std::cout << "x: " << p.x << std::endl;
  std::cout << "y: " << p.y << std::endl;
  std::cout << "theta: " << p.theta << std::endl;
  std::cout << "weight: " << p.weight << " }" << std::endl;
  std::cout << std::endl;
}

void ParticleFilter::PrintLandmarkObs(const LandmarkObs& obs) const
{
  std::cout << "LandmarkObs: { " << std::endl;
  std::cout << "id: " << obs.id << std::endl;
  std::cout << "x: " << obs.x << std::endl;
  std::cout << "y: " << obs.y << std::endl;
  std::cout << std::endl;
}
