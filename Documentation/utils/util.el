
(defun edit-mlx (mlxfilename editfunc)
  "Unpack mlxfilename in a the directory Temp (should be created
apriori), run the function editfunc in the matlab/document.xml
file, pack the result again and delete the unpacked files in
Temp."
  (let ((command-str ""))
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
    (setq command-str (concat command-str "mv temp.mlx .. \n"))
    (setq command-str (concat command-str "cd .. \n"))
    (setq command-str (concat command-str "rm -r Temp/* \n"))
    (shell-command command-str)
    )
  )


;; example
(defun editfunc ()
  (while (re-search-forward (rx "Introduction") nil t)
    (replace-match "hello")
    )
  )

;; example
;; (edit-and-compress "tutorial_2_changing_control_protocol_live.mlx" #'editfunc)


