import { Component, OnInit, Input } from '@angular/core';
import { ModalController } from '@ionic/angular';
import { SafeResourceUrl } from '@angular/platform-browser';

@Component({
  selector: 'app-show-mirna-modal',
  templateUrl: './show-mirna-modal.component.html',
  styleUrls: ['./show-mirna-modal.component.scss'],
})
export class ShowMirnaModalComponent implements OnInit {
  @Input() mirna: string;
  @Input() document: SafeResourceUrl;

  constructor(private modalCtrl: ModalController) { }

  ngOnInit() { }

  onCloseModal() {
    this.modalCtrl.dismiss();
  }

}
