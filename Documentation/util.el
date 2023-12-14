;;; Utils for building BattMo doc

;; Use user-login-name to get user name
;; (cond '(compare-strings user-login-name "xavier")

(pcase (user-login-name)
  ("xavier" (progn
              (pyvenv-activate "~/Python/battmodoc-env/")
              (setq docdir "/home/xavier/Matlab/Projects/battmo/Documentation/")
              (setq utilsdir "/home/xavier/Matlab/Projects/battmo/Documentation/utils/")
              (setq testdir "/home/xavier/Matlab/Projects/battmo-doc-test/")))
  ("august" (progn
              (setq docdir "/home/august/Projects/Battery/BattMo-dev/Documentation/")
              (setq utilsdir "/home/august/Projects/Battery/BattMo-dev/Documentation/utils/")
              (setq testdir "/home/august/Projects/Battery/battmo-doc-test/")
              
              )
   )
  )


     
(defun battmodoc-build ()
  (interactive)
  (let ((outputbuffer (get-buffer-create "*buildoutput*"))
        (default-directory docdir))
    (pop-to-buffer outputbuffer)
    (erase-buffer)
    (start-process "battmo-build" outputbuffer "make" "html")
    )
  (browse-url (concat docdir "_build/html/index.html"))
  )

(defun battmodoc-build-examples ()
  (interactive)
  (let ((outputbuffer (get-buffer-create "*publishoutput*"))
        (directory utilsdir))
    (pop-to-buffer outputbuffer)
    (erase-buffer)
    (start-process "battmo-publish" outputbuffer "python" (concat directory "buildPublishedExamples.py"))
    )
  )

(defun battmodoc-publish ()
  (interactive)
  (let ((default-directory docdir))
    (shell-command (concat "cp -rf _build/html/* " testdir))
    )
  (let ((default-directory testdir))
    (shell-command "git add *" )
    (shell-command "git commit -m \"update in doc\"")
    (shell-command "git push -f" )
    )
  (browse-url "https://github.com/BattMoTeam/battmo-doc-test/actions")
  (browse-url "https://battmoteam.github.io/battmo-doc-test/")
  )

(defun convert-to-attribute ()
  (interactive)
  (re-search-forward (rx (group (1+ (any word))) (0+ space) "%" (0+ space) (group (0+ nonl))))
  (let ((att (match-string 1))
        (com (match-string 2))
        )
    (beginning-of-line)
    (kill-line)
    (insert "   .. attribute:: " att "\n\n")
    (insert "      " com "\n")
    )
  )

