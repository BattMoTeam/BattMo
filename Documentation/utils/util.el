
(defun edit-mlx (mlxfilename editfunc)
  "Unpack mlxfilename in a the directory Temp (should be created
apriori), run the function editfunc in the matlab/document.xml
file, pack the result again and delete the unpacked files in
Temp."
  (let ((mlxfilename (expand-file-name mlxfilename))
        (command-str ""))
    (setq command-str (concat "unzip " mlxfilename " -d Temp\n"))
    (shell-command command-str)
    (with-temp-buffer
      (insert-file "Temp/matlab/document.xml")
      (beginning-of-buffer)
      (funcall editfunc)
      (write-region (point-min) (point-max) "Temp/matlab/document.xml")
      )
    (setq command-str "cd Temp \n")
    (setq command-str (concat command-str "zip -r temp.mlx * \n"))
    (setq command-str (concat command-str "mv -f temp.mlx " mlxfilename " \n"))
    (setq command-str (concat command-str "cd .. \n"))
    (setq command-str (concat command-str "rm -r Temp/* \n"))
    (shell-command command-str)
    )
  )


;; example
(defun editfunc ()
  (while (re-search-forward (rx "plotBatteryMesh") nil t)
    (replace-match "plotBatteryGrid")
    )
  )

;; example
(edit-mlx "../../Examples/Notebooks/tutorial_7_a_simple_p4d_model_live.mlx" #'editfunc)


