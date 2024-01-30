;;; Utils for building BattMo doc

;; Use user-login-name to get user name
;; (cond '(compare-strings user-login-name "xavier")

(pcase (user-login-name)
  ("xavier" (progn
              (pyvenv-activate "~/Python/battmodoc-env/")
              (setq docdir "/home/xavier/Matlab/Projects/battmo/Documentation/")
              (setq rootdir "/home/xavier/Matlab/Projects/battmo/")
              (setq utilsdir "/home/xavier/Matlab/Projects/battmo/Documentation/utils/")
              (setq testdir "/home/xavier/Matlab/Projects/battmo-doc-test/")))
  ("august" (progn
              (setq docdir "/home/august/Projects/Battery/BattMo/Documentation/")
              (setq utilsdir "/home/august/Projects/Battery/BattMo/Documentation/utils/")
              (setq testdir "/home/august/Projects/Battery/battmo-doc-test/")

              )
   )
  )

(defun battmodoc-local-open ()
  "Open locally built documentation in browser"
  (interactive)
  (browse-url (concat docdir "_build/html/index.html"))
  )

(defun battmodoc-build ()
  "Build BattMo documentation"
  (interactive)
  (let* ((default-directory docdir)
         (outputbuffer (get-buffer-create "*buildoutput*"))
         )
    (pop-to-buffer outputbuffer)
    (erase-buffer)
    (start-process "battmo-build" outputbuffer "make" "html")
    )
  )

(defun battmodoc-build-examples ()
  "Run python publish script"
  (interactive)
  (let ((outputbuffer (get-buffer-create "*publishoutput*"))
        (directory utilsdir))
    (pop-to-buffer outputbuffer)
    (erase-buffer)
    (start-process "battmo-publish" outputbuffer "python" (concat directory "buildPublishedExamples.py"))
    )
  )

(defun battmodoc-publish ()
  "Copy and paste build directory in doc repo and push the result. Then open in a browser the action page in the repo and the published page"
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

(defun battmodoc-matlab-include-link ()
  (interactive)
  (let* ((url "https://github.com/BattMoTeam/BattMo")
        (branch-name "dev")
        (fullfilename (battmo-counsel-fzf "" rootdir))
        (filename (file-name-nondirectory fullfilename))
        (p (make-marker)))
    (insert (concat "<" url "/blob/" branch-name "/" fullfilename "  "))
    (set-marker p (point))
    (insert (concat filename ">"))
    (goto-char p)
    )
  )

(defun battmo-counsel-fzf (initial-input initial-directory)
  (counsel-require-program counsel-fzf-cmd)
  (setq counsel--fzf-dir  initial-directory)
  (ivy-read "file to link:"
            #'counsel-fzf-function
            :initial-input initial-input
            :re-builder #'ivy--regex-fuzzy
            :dynamic-collection t
            :action #'identity
            :caller 'counsel-fzf))

