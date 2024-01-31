(defun capture-field ()
  (interactive)
  (let ((name (thing-at-point 'word))
        (str (save-excursion
               (re-search-forward "{" nil t)
               (let* ((start (progn (backward-char)
                                    (point)))
                      (end (progn (xah-goto-matching-bracket)
                                  (point)
                                  )))
                 (buffer-substring-no-properties start end)
                 ))))
    (let-alist (json-read-from-string str)
      (kill-new (concat name " % " .description " (symbol: " .symbol ")"))
      )
    )
  )

(let-alist (with-current-buffer "SolidDiffusionModel.schema.json"
             (beginning-of-buffer)
             (json-parse-buffer :object-type 'alist))
  .properties.particleRadius.description)




(json-read)
