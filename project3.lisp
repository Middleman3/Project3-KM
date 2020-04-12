;; A K-Means implementation

;; some handy parameters
(defparameter *debug-level* 2 "how many debugs do you wish?")
(defparameter *p* 2 "Used as exponent in distance formula")
(defparameter *point-dimensions* 2 "The dimension of the data points")
(defparameter *max-value* 1000 "The maximum value of a data point component")
(defparameter *stopping-friction* 5 "Sum centroid motion needed to not be considered converged")

;;;;;;;;;;;;;;; Some Recycled Functions

(defun display (debug-level description &rest args)
  (if (<= debug-level *debug-level*)
      (apply #'format t description args)))

(defmacro while (test return-value &body body)
  "Repeatedly executes body as long as test returns true.
Then evaluates and returns return-value"
  `(progn (loop while ,test do (progn ,@body)) ,return-value))

(defun random-elt (sequence)
  "Returns a random element from sequence"
  (elt sequence (random (length sequence))))

(defun generate-list (num function &optional no-duplicates)
  "Generates a list of size NUM, with each element created by
  (funcall FUNCTION).  If no-duplicates is t, then no duplicates
are permitted (FUNCTION is repeatedly called until a unique
new slot is created).  EQUALP is the test used for duplicates."
  (let (bag)
    (while (< (length bag) num) bag
      (let ((candidate (funcall function)))
	(unless (and no-duplicates
		     (member candidate bag :test #'equalp))
	  (push candidate bag))))))

;;;;;;;; New Generic Helper Functions

(defun position-of-best (seq fitness aggregate)
  "Returns the best index of the seq in regards to the given fitness function and aggregation method (min, max, perhaps medoid)"
  (let* ((fitnesses (mapcar fitness seq)))
    (position (apply aggregate fitnesses) fitnesses)))

(defun default (vals defaults)
  "given a parallel list of values and default values, returns a list of values, where nils defaulted to the default value"
  (mapcar (lambda (val def) (or val def)) vals defaults))

;;;;;;;; K Means Sub Functions

(defun lp-distance (point1 point2 &optional (p *p*))
  "Calculates the l^p distance between two n-dimensional points. p=1 -> manhattan, p=2 -> Euclidean, etc"
  (expt (apply #'+ (mapcar (lambda (e1 e2)
			     (expt (abs (- e1 e2)) p))
			   point1
			   point2))
	(/ 1 p)))	    

(defun my-centroids-index (point centroids)
  "Uses lp-distance to find the closest centroid"
  (position-of-best centroids (lambda (centriod) (lp-distance point centriod)) #'min))

(defun mean-of-cluster (cluster)
  "Calculates the mean of a cluster of points. Assumes points follow *point-dimensions*"
  (if (not cluster) nil ; if cluster is empty, return nil
      (let ((sum-vector (make-list *point-dimensions* :initial-element 0))
	    (size (length cluster)))
	(mapcar (lambda (point) (setf sum-vector (mapcar #'+ point sum-vector))) cluster) ; sum all points into sum-vector
	(mapcar (lambda (sum) (/ sum size)) sum-vector)))) ; divide sum-vector by size of cluster and return    

(defun transpose (matrix)
  "Transposes a list of lists matrix. Helps with find-disparate-points"
  (apply #'mapcar #'list matrix))

#|
(defun find-disparate-points (points num-clusters)
  "Analyzes points to select distant initial clusters. Checks the min, max, and closest to zero of each dimension among points."
  (let* ((count 0)
	 (methods (list (lambda (x) x) ; to find min value
			(lambda (x) (- x)) ; to find max value
			(lambda (x) (abs x)))) ; to find value closest to 0
	 (hits (* (length methods) *point-dimensions*)) ; hits = how many values I can come up with that are disparate
	 (contests (transpose points)) ; group points by components: ((x1 y1) (x2 y2)) -> ((x1 x2) (y1 y2))
	 index ; index of the point that won the current constest
	 result)
    (labels ((choose ()) ;;; PICK UP HERE. 
    (generate-list num-clusters
		   (lambda ()
		     (if (< i hits) ; Unsure how many centroids I'll need to initialize, but I have a plan for the first 'hits centriods
			 (multiple-value-bind (j k) (floor i (length methods)) ; i gives me the k'th method on j'th contest
			   (setf index (position-of-best (elt contests j) (elt methods k) #'min)) ; save index of winner
			   (setf result (append result (list (elt points index))))) ; append winning point to result because its disparate
			 (setf result (append result (list (random-elt points)))))))) ; out of ideas for disparity and still need more centroids)
    
|#

(defun initialize-centroids (points num-clusters)
  "Returns a new set of centroids. Currently chooses random points."
  (generate-list num-clusters (lambda () (random-elt points)) t))

(defun new-centroids (clusters)
  "Given a list of clusters, returns a list of the new centroids"
  (let (centroids)
    (dolist (cluster clusters centroids)
      (setf centroids (append centroids (list (mean-of-cluster cluster))))))) ; append new centroid to centroids list

(defun assign-points (points centroids)
  "Given a set of points and a list of centroids, returns a list of clusters, parallel to the list of centroids"
  (let ((clusters (make-list (length centroids)))
	index)
    (dolist (point points clusters)
      (setf index (my-centroids-index point centroids)) ; find closest centroid
      (setf (elt clusters index) (append (elt clusters index) (list point)))))) ; append point to closest centroid's cluster

(defparameter *motion* 0 "Sum measure of change in centroids")

(defun continue-p (prev-centroids new-centroids)
  "Compare the motion of each centroid to determine if they have collectively converged. Returns t if execution should continue"
  (setf *motion* (apply #'+ (mapcar #'lp-distance prev-centroids new-centroids)))
  (> *motion* *stopping-friction*))
       
(defun csv-write-points (points)
  "comma separates a list of lists"
  (dolist (point points)
    (dolist (val point) (format t "~D," val))
    (format t "~%")))

(defun csv-write-centroids (centroids iteration)
  "writes csv to be appended to points, and plotted in Google Sheets"
  (dolist (centroid centroids)
    (format t "~F," (first centroid))
    (dotimes (i iteration) (format t ","))       
    (dolist (val (rest centroid)) (format t "~F," val))
    (format t "~%")))

(defun k-means-clustering (points num-clusters &optional to-csv)
  "Given a set of data points and a number of clusters, performs k-means clustering, returning a list of cluster centroids."
  (if to-csv (csv-write-points points))
  (let ((centroids (initialize-centroids points num-clusters))       
	prev-centroids
	clusters
	(count 0))
    (while (or (continue-p prev-centroids centroids) (= count 0)) centroids
      (if to-csv	  
	  (csv-write-centroids centroids (incf count))
	  (display 0 "~%------ Iteration ~D | motion = ~D -----~%" (incf count) *motion*))
      (setf prev-centroids centroids) ; shift perspective of time
      (setf clusters (assign-points points centroids)) ; determine clusters from centroids
      (setf centroids (default (new-centroids clusters) centroids))))) ; Alter centroids, don't change centroids of empty clusters     

;;;;;;;;;;;;;;;; dont forget to change back to gaussian random
