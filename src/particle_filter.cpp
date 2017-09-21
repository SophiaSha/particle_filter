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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;

	num_particles = 20;

	particles.clear();
	weights.clear();

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	for(int i=0; i < num_particles ;i++){

		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1/num_particles;

		particles.push_back(particle);
		weights.push_back(particle.weight);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	for(int i=0; i < num_particles ;i++){
		if(fabs(yaw_rate)<0.0001){
			particles[i].x = particles[i].x + velocity*cos(particles[i].theta)*delta_t;
			particles[i].y = particles[i].y + velocity*sin(particles[i].theta)*delta_t;
		}
		else{
			particles[i].x = particles[i].x + (sin(particles[i].theta+delta_t*yaw_rate)-sin(particles[i].theta))*velocity/yaw_rate;
			particles[i].y = particles[i].y + (cos(particles[i].theta)-cos(particles[i].theta+delta_t*yaw_rate))*velocity/yaw_rate;
			particles[i].theta = particles[i].theta+delta_t*yaw_rate;
		}

	normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
	normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
	normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

	particles[i].x = dist_x(gen);
	particles[i].y = dist_y(gen);
	particles[i].theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	double min_dist;
	double temp_dist;
	LandmarkObs obs;
	LandmarkObs pred;
	double target_x;
	double target_y;
	int ID=-1;

	for(int i=0; i < observations.size() ;i++){
		obs = observations[i];
		min_dist = numeric_limits<double>::max();
		for(int j=0; j < predicted.size();j++){
			pred = predicted[j];
			temp_dist = dist(obs.x, obs.y, pred.x, pred.y);
			if(temp_dist<min_dist){
				min_dist = temp_dist;
				ID = j;
			}
		observations[i].x = predicted[ID].x;
		observations[i].y = predicted[ID].y;
		}
	}
}

double WeightCalcu(double sig_x, double sig_y, double x_obs, double y_obs, double mu_x, double mu_y){

	// calculate normalization term
	//double gauss_norm= (1/(2 * M_PI * sig_x * sig_y));

	// calculate exponent
	//double exponent= pow((x_obs - mu_x),2)/(2 * pow(sig_x,2)) + pow((y_obs - mu_y),2)/(2 * pow(sig_y,2));

	// calculate weight using normalization terms and exponent
	double weight= (1/(2 * M_PI * sig_x * sig_y)) * exp(-(pow((x_obs - mu_x),2)/(2 * pow(sig_x,2)) + pow((y_obs - mu_y),2)/(2 * pow(sig_y,2))));

	return weight;
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
	double x;
	double y;
	double theta;
	//cout<<"map.size:" <<map_landmarks.landmark_list.size()<<endl;
	for(int i=0;i < num_particles;i++){
		 x = particles[i].x;
		 y = particles[i].y;
		 theta = particles[i].theta;

		 vector<LandmarkObs> predicted; //predicted in MAP'S coordinate system
		 vector<LandmarkObs> obs_map;   //observations in MAP'S coordinate system

		 for(int m=0;m<map_landmarks.landmark_list.size();m++){
			 if(dist(x,y,map_landmarks.landmark_list[m].x_f,map_landmarks.landmark_list[m].y_f) <= sensor_range*sensor_range){
				 LandmarkObs candi;
				 candi.id = map_landmarks.landmark_list[m].id_i;
				 candi.x = map_landmarks.landmark_list[m].x_f;
				 candi.y = map_landmarks.landmark_list[m].y_f;
				 predicted.push_back(candi);
			 }
		 }
		 //cout<<"observations.size:"<<observations.size()<<endl;

		 for(int m=0;m<observations.size();m++){
			 	LandmarkObs obs;
				obs.id= observations[m].id;
				obs.x= x + (cos(theta) * observations[m].x) - (sin(theta) * observations[m].y);
				obs.y= y + (sin(theta) * observations[m].x) + (cos(theta) * observations[m].y);
				obs_map.push_back(obs);
		 }
		 vector<LandmarkObs> obs_map_backup = obs_map;
		 dataAssociation(predicted,obs_map);

		 particles[i].weight = 1.0;
	 
		 //cout<<"predicted.size:"<<endl;
		 for(int m=0;m<obs_map_backup.size();m++){
			 double w = WeightCalcu(std_landmark[0], std_landmark[1], obs_map_backup[m].x, obs_map_backup[m].y, obs_map[m].x, obs_map[m].y);
			 particles[i].weight = particles[i].weight * w;
			 
		 }
		 weights[i] = particles[i].weight;
		 
	}


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	default_random_engine gen;
	vector<Particle> resample_particles;

	discrete_distribution<> d(weights.begin(),weights.end());
	
	for(int i=0;i<particles.size();i++){
		resample_particles.push_back(particles[d(gen)]);
	}
	particles = resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
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
