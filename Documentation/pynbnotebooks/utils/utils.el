(setq filenames (list "tutorial_1_a_simple_p2d_model_live"
                      "tutorial_2_changing_control_protocol_live"
                      "tutorial_3_modify_structural_parameters_live"
                      "tutorial_4_modify_material_parameters_live"
                      "tutorial_5_simulate_CCCV_cycling_live"
                      "tutorial_6_simulate_thermal_performance_live"
                      "tutorial_7_a_simple_p4d_model_live"
                      "tutorial_8_simulate_a_multilayer_pouch_cell_live"
                      "tutorial_9_simulate_a_cylindrical_cell_live"))


(defun clean-up-ipynb-output ()
  "removes the verbose output from the simulations from a ipynb notebook"
  (interactive)
  (let* ((find-next-output-region (lambda ()
                                   (if (re-search-forward "outputs\": " nil t)
                                       (save-excursion
                                         (let ((start (point))
                                               (end (progn (forward-sexp)
                                                           (point))))
                                           (list start end)
                                           ))
                                     nil
                                     )))
         (res (funcall find-next-output-region))
         )
    (while res
      (if (and (re-search-forward "data" nil t)
               (re-search-forward "Solving timestep" nil t)
               (< (point) (cadr res)))
          (progn (goto-char (car res))
                 (delete-region (car res) (cadr res))
                 (insert "[]"))
        (goto-char (cadr res)))
      (setq res (funcall find-next-output-region)))
    )
  (save-buffer)
  )

(defun get-details (filename)
  (with-current-buffer (find-file-noselect (concat "/home/xavier/Matlab/Projects/battmo/Documentation/utils/Temp/" filename ".m"))
    (beginning-of-buffer)
    (re-search-forward (rx "Tuto" (0+ (not numeric)) (group (any numeric)) (0+ (not "-")) "- " (group (1+ any)) eol))
    `((filename . ,filename)
      (num . ,(match-string-no-properties 1))
      (title . ,(match-string-no-properties 2)))))

(defun setup-page (notebook)
  (with-temp-buffer 
    (erase-buffer)
    (let-alist notebook
      (insert ".. _tutorial" .num ":\n\n" )
      (insert .title "\n")
      (insert (make-string (length .title) ?=) "\n\n")
      (insert "The notebook is available in your BattMo installation. Run\n\n")
      (insert ".. code:: matlab\n\n")
      (insert "   open " .filename "\n\n")
      (insert ".. raw:: html\n")
      (insert "   :file: ../_static/notebooks/" .filename ".html\n")
      )
    (buffer-string)
    )
  )

(defun create-page (filename)
  (let* ((directory "/home/xavier/Matlab/Projects/battmo/Documentation/tutorials/")
         (notebook (get-details filename))
         (content (setup-page notebook))
         (outputname (concat directory "tutorial" (let-alist notebook .num) ".rst")))
    (with-current-buffer (find-file-noselect outputname)
      (erase-buffer)
      (insert content)
      (save-buffer)
      )
    )
  )

(cl-loop for filename in filenames do (create-page filename))

(let ((content (cl-loop for filename in filenames concat
                        (let-alist (get-details filename)
                          (concat "   " .title " <tutorials/tutorial" .num ">\n")))))
  (with-current-buffer "tutorials.rst"
    (beginning-of-buffer)
    (re-search-forward ":hidden:")
    (insert "\n\n")
    (insert content)
    ))


(let ((content (cl-loop for filename in filenames concat
                        (let-alist (get-details filename)
                          (concat "   .. grid-item-card::\n"
                                  "      :padding: 2\n"
                                  "\n"
                                  "      :ref:`" .title "<" filename ".pynb>`\n\n")))))
  (with-current-buffer "tutorials.rst"
    (beginning-of-buffer)
    (re-search-forward ".. grid:: 1")
    (insert "\n\n")
    (insert content)
    ))


(let ((content (cl-loop for filename in filenames concat
                        (let-alist (get-details filename)
                          (concat "   .. grid-item-card::\n"
                                  "      :padding: 2\n"
                                  "\n"
                                  "      :ref:`" .title "<tutorial" .num ">`\n\n")))))
  (with-current-buffer "tutorials.rst"
    (beginning-of-buffer)
    (re-search-forward ".. grid:: 1")
    (insert "\n\n")
    (insert content)
    ))

