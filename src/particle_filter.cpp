/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

#define EPS 0.0001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles=100;

	//std::vector<double> weights(num_particles, 1.0f);

	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {
		Particle p;
		p.id=i;
		p.x +=dist_x(gen);
		p.y += dist_y(gen);
		p.theta += dist_theta(gen);
		p.weight=1.0;
		particles.push_back(p);
	}

	is_initialized=true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	for (int i = 0; i < num_particles; ++i) {
		double theta=particles[i].theta;
		double x=particles[i].x;
		double y=particles[i].y;

		double x_new=x+velocity*delta_t*cos(theta);
		double y_new = y+velocity*delta_t*sin(theta);
		double theta_new = theta;
		if (fabs(yaw_rate) > EPS) {
			x_new=x+velocity*(sin(theta+delta_t*yaw_rate)-sin(theta))/yaw_rate;
			y_new = y+velocity*(cos(theta)-cos(theta+delta_t*yaw_rate))/yaw_rate;
			theta_new = theta+delta_t*yaw_rate;
		}
		normal_distribution<double> dist_x(x_new, std_pos[0]);
		normal_distribution<double> dist_y(y_new, std_pos[1]);
		normal_distribution<double> dist_theta(theta_new, std_pos[2]);

		particles[i].x=dist_x(gen);
		particles[i].y=dist_y(gen);
		particles[i].theta=dist_theta(gen);
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for (unsigned int i = 0; i < observations.size(); i++) { //look for measurement closest to landmark

	    double dist_min = numeric_limits<double>::max(); // init closest distance
	    int idx = -1; //init index of closest measurement

	    for (unsigned int j = 0; j < predicted.size(); j++) {
	      LandmarkObs measurement = predicted[j];
	      double dist_current = dist(observations[i].x, observations[i].y, measurement.x, measurement.y);

	      // find the predicted landmark nearest the current observed landmark
	      if (dist_current < dist_min) {
	    	  dist_min = dist_current;
	    	  idx = measurement.id;
	      }
	    }
	    observations[i].id = idx;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for(int i = 0; i < num_particles; ++i){
		double theta=particles[i].theta;
		double x=particles[i].x;
		double y=particles[i].y;
		///1st step: transform sensor measurements to the map coordinate system (absolute system)
	    vector<LandmarkObs> observations_abs;
	    for (unsigned int j = 0; j < observations.size(); j++) {
	      double x_abs = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
	      double y_abs = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
	      observations_abs.push_back(LandmarkObs{observations[j].id, x_abs, y_abs });
	    }
	    ///2nd step: keep landmarks within sensor_range distance of the particle
	    vector<LandmarkObs> inRange;
	    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
	      double x_l = map_landmarks.landmark_list[j].x_f;
	      double y_l = map_landmarks.landmark_list[j].y_f;
	      int id_l = map_landmarks.landmark_list[j].id_i;
	      double dist_landmark=dist(x,y,x_l,y_l);
	      if (dist_landmark<= sensor_range) {
	    	  inRange.push_back(LandmarkObs{id_l, x_l, y_l});
	      }
	    }
		///3nd step: associate measurements to the landmarks "seen" by particle
	    dataAssociation(inRange, observations_abs);
		///4th step: update weight
	    particles[i].weight = 1.0;
	    double s_x = std_landmark[0];
	    double s_y = std_landmark[1];
	    //go through each measurement "taken" by the particle
	    for (unsigned int j = 0; j < observations_abs.size(); j++) {
	      double meas_x = observations_abs[j].x;
	      double meas_y = observations_abs[j].y;
	      double x_l,y_l;
	      // find closest landmark coordinates
	      for (unsigned int k = 0; k < inRange.size(); k++) {
	        if (inRange[k].id == observations_abs[j].id) {
	          x_l = inRange[k].x;
	          y_l = inRange[k].y;
	        }
	      }
	      //calculate measurement weight
	      double norm=1/(2*M_PI*s_x*s_y);
	      double gaussian_term1=pow(meas_x-x_l,2)/(2*pow(s_x, 2));
	      double gaussian_term2=pow(meas_y-y_l,2)/(2*pow(s_y, 2));
	      double meas_w =norm*exp(-gaussian_term1-gaussian_term2);
	      //multiply it to other measurement's weights
	      particles[i].weight *= meas_w;
	    }
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	vector<Particle> new_particles;
	vector<double> weights;
	default_random_engine gen;
	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight);
	}
	double max_weight = *max_element(weights.begin(), weights.end());

	uniform_real_distribution<double> dist_reel(0.0, max_weight);
	uniform_int_distribution<int> dist_discrete(0, num_particles - 1);

	int idx = dist_discrete(gen);
	double beta = 0.0;

	for (int i = 0; i < num_particles; i++) {
		beta += dist_reel(gen) * 2.0;
	    while (beta > weights[idx]) {
	      beta -= weights[idx];
	      idx = (idx + 1) % num_particles;
	    }
	    new_particles.push_back(particles[idx]);
	}
	particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
