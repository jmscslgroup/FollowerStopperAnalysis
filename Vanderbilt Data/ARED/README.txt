Arizona Ring Experiments Dataset

The Arizona Ring Experiments Dataset (ARED) holds the measurements for the
eight experiments described in the following article:

F. Wu, R. Stern, S. Cui, M. L. Delle Monache, R. Bhadani, M. Bunting,
M. Churchill, N. Hamilton, R. Haulcy, B. Piccoli, B. Seibold, J. Sprinkle,
D. Work. “Tracking vehicle trajectories and fuel rates in oscillatory traffic.”
to appear in Transportation Research Part C: Emerging Technologies. 2017.

If our dataset is helpful for your research, please consider citing our work
For inquiry or to report a bug, please contact the author Fangyu Wu
through fangyuwu(at)berkeley(dot)edu.

## Specifications
The definition of each data attribute is explained below:
-exp_id: experiment identification number
-veh_id: vehicle identification number
-length: vehicle length (m)
-height: vehicle height (m)
-make: vehicle make
-year: vehicle year
-model: vehicle model
-timestamp: timestamp of this measurement (s)
-cam_distance: distance traveled from the origin of the camera frame (m)
-cam_velocity: vehicle velocity estimated from the camera (m/s)
-cam_acceleration: vehicle acceleration estimated from the camera (m/s2)
-leader_id: leading vehicle identification number
-cam_leader_gap: space headway from current vehicle to leading vehicle
estimated from the camera (m)
-follower_id: following vehicle identification number
-cam_follower_gap: space headway from current vehicle to following vehicle
estimated from the camera (m)
-obd_velocity: vehicle velocity measured by OBD (m/s)
-obd_fuel_rate: vehicle fuel rate measured by OBD (l/h)
-cam_upper_front_x: _unsmoothed_ x coordinate of the upper front corner of
vehicle's bounding box estimated from the camera (m)
-cam_upper_front_y: _unsmoothed_ y coordinate of the upper front corner of
vehicle's bounding box estimated from the camera (m)
-cam_lower_rear_x: _unsmoothed_ x coordinate of the lower rear corner of
vehicle's bounding box estimated from the camera (m)
-cam_lower_rear_y: _unsmoothed_ y coordinate of the lower rear corner of
vehicle's bounding box estimated from the camera (m)

Note that "nan" entries indicate that either the data is missing or is invalid.
For more information, please refer to our
[manuscript](https://hal.archives-ouvertes.fr/hal-01614665).

The work be cited as follows:

@article{wu2019tracking,
  title={Tracking vehicle trajectories and fuel rates in oscillatory traffic},
  author={Wu, F. and Stern, R. and Cui, S. and Delle Monache, M. L. and Bhadani, R. and Bunting, M. and Churchill, M. and Hamilton, N. and Haulcy, R. and Pohlmann, H. and Piccoli, B. and Seibold, B. and Sprinkle, J. and Work, D. B.},
  journal={Transportation Research Part C: Emerging Technologies, to appear},
  year={2019},
  publisher={Elsevier}
}

## License

Arizona Ring Experiments Dataset (c) by Work Research Group

Arizona Ring Experiments Dataset is licensed under a
Creative Commons Attribution-ShareAlike 3.0 Unported License.
You should have received a copy of the license along with this
work.  If not, see <http://creativecommons.org/licenses/by-sa/3.0/>.

Last updated on 01/04/2019